#! /usr/bin/env python2.7
# -*- coding: utf-8 -*-
#
# pySpatialETL is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# any later version.
#
# pySpatialETL is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# Author : Fabien Rétif - fabien.retif@zoho.com
#
from __future__ import division, print_function, absolute_import

import os
from itertools import product
from multiprocessing import Pool, cpu_count

import geopandas as gpd
import numpy as np
import pandas as pd
from osgeo import gdal, osr, gdal_array
from osgeo import ogr
from rasterio.crs import CRS
from shapely.ops import unary_union

from spatialetl.coverage.TimeCoverage import TimeCoverage
from spatialetl.coverage.TimeLevelCoverage import TimeLevelCoverage
from spatialetl.coverage.io.CoverageWriter import CoverageWriter
from spatialetl.exception.CoverageError import CoverageError
from spatialetl.utils.VariableDefinition import VariableDefinition
from spatialetl.utils.logger import logging
from sklearn.cluster import AgglomerativeClustering

import numpy as np
from shapely.ops import split
import geopandas
from shapely.geometry import MultiPolygon, Polygon, LineString

def union(x):
    return unary_union(x)




def rhombus(square):
    """
    Naively transform the square into a Rhombus at a 45 degree angle
    """
    coords = square.boundary.coords.xy
    xx = list(coords[0])
    yy = list(coords[1])
    radians = 1
    points = list(zip(xx, yy))
    Rhombus = Polygon(
        [
            points[0],
            points[1],
            points[3],
            ((2 * points[3][0]) - points[2][0], (2 * points[3][1]) - points[2][1]),
            points[4],
        ]
    )
    return Rhombus


def get_squares_from_rect(RectangularPolygon, side_length=0.0025):
    """
    Divide a Rectangle (Shapely Polygon) into squares of equal area.

    `side_length` : required side of square

    """
    rect_coords = np.array(RectangularPolygon.boundary.coords.xy)
    y_list = rect_coords[1]
    x_list = rect_coords[0]
    y1 = min(y_list)
    y2 = max(y_list)
    x1 = min(x_list)
    x2 = max(x_list)
    width = x2 - x1
    height = y2 - y1

    xcells = int(np.round(width / side_length))
    ycells = int(np.round(height / side_length))

    yindices = np.linspace(y1, y2, ycells + 1)
    xindices = np.linspace(x1, x2, xcells + 1)
    horizontal_splitters = [
        LineString([(x, yindices[0]), (x, yindices[-1])]) for x in xindices
    ]
    vertical_splitters = [
        LineString([(xindices[0], y), (xindices[-1], y)]) for y in yindices
    ]
    result = RectangularPolygon
    for splitter in vertical_splitters:
        result = MultiPolygon(split(result, splitter))
    for splitter in horizontal_splitters:
        result = MultiPolygon(split(result, splitter))
    square_polygons = list(result)

    return square_polygons


def split_polygon(G, side_length=0.025, shape="square", thresh=0.9):
    """
    Using a rectangular envelope around `G`, creates a mesh of squares of required length.

    Removes non-intersecting polygons.


    Args:

    - `thresh` : Range - [0,1]

        This controls - the number of smaller polygons at the boundaries.

        A thresh == 1 will only create (or retain) smaller polygons that are
        completely enclosed (area of intersection=area of smaller polygon)
        by the original Geometry - `G`.

        A thresh == 0 will create (or retain) smaller polygons that
        have a non-zero intersection (area of intersection>0) with the
        original geometry - `G`

    - `side_length` : Range - (0,infinity)
        side_length must be such that the resultant geometries are smaller
        than the original geometry - `G`, for a useful result.

        side_length should be >0 (non-zero positive)

    - `shape` : {square/rhombus}
        Desired shape of subset geometries.


    """
    assert side_length > 0, "side_length must be a float>0"
    Rectangle = G.envelope
    squares = get_squares_from_rect(Rectangle, side_length=side_length)
    SquareGeoDF = geopandas.GeoDataFrame(squares).rename(columns={0: "geometry"})
    Geoms = SquareGeoDF[SquareGeoDF.intersects(G)].geometry.values
    if shape == "rhombus":
        Geoms = [rhombus(g) for g in Geoms]
        geoms = [g for g in Geoms if ((g.intersection(G)).area / g.area) >= thresh]
    elif shape == "square":
        geoms = [g for g in Geoms if ((g.intersection(G)).area / g.area) >= thresh]
    return geoms


class SnaleWriter (CoverageWriter):

    def __init__(self,cov,myFile):
        CoverageWriter.__init__(self,cov,myFile);

        if self.coverage.is_regular_grid() == False:
            raise ValueError("This writer supports only Coverage with a regular horizontal axis.")

        if os.path.isdir(self.filename) is False:
            raise ValueError("Filename has to be a directory.")

        gdal.AllRegister()
        self.driver = ogr.GetDriverByName("ESRI Shapefile")

        self.rows = self.coverage.get_x_size(type="target_global")
        self.cols = self.coverage.get_y_size(type="target_global")

        xmin=np.min(self.coverage.read_axis_x(type="target_global"))
        xmax=np.max(self.coverage.read_axis_x(type="target_global"))
        ymin=np.min(self.coverage.read_axis_y(type="target_global"))
        ymax=np.max(self.coverage.read_axis_y(type="target_global"))
        self.x_pixel_size=round((xmax-xmin)/self.coverage.get_x_size(type="target_global"),6)
        self.y_pixel_size=round((ymax-ymin)/self.coverage.get_y_size(type="target_global"),6)

        self.geotransform = (xmin, self.x_pixel_size, 0, ymin, 0, self.y_pixel_size)

    def close(self):
        return

    # Variables

    #################
    # HYDRO
    # Sea Surface
    #################
    def write_variable_sea_water_column_thickness(self):

        if (isinstance(self.coverage, TimeCoverage) or isinstance(self.coverage, TimeLevelCoverage)):

            # if self.coverage.rank == 0:
            #     logging.info('[DefaultWriter] Gathering data from all processors')
            #     global_data = np.empty(
            #         [self.coverage.get_t_size(type="target_global"),
            #          self.coverage.get_y_size(type="target_global"),
            #          self.coverage.get_x_size(type="target_global")])
            #     global_data[:] = np.nan

            for time_index in range(0, self.coverage.get_t_size()):
                time = self.coverage.read_axis_t(type="target_global")[time_index]
                local_data = self.coverage.read_variable_sea_water_column_thickness_at_time(time_index)

                local_data[local_data == 9.96921e+36] = 0  # Remove NaN
                local_data[local_data <= 0.1] = 0
                local_data[(local_data > 1.2) & (local_data != 0)] = 3
                local_data[(local_data >= 0.7) & (local_data <= 1.2)] = 2
                local_data[(local_data > 0.1) & (local_data <= 0.7)] = 1
                #local_data[(local_data > 0.1) & (local_data <= 0.3)] = 1

                x, y, elevation = self.coverage.read_axis_x(type="target_local"), self.coverage.read_axis_y(
                    type="target_local"), local_data

                x, y = np.meshgrid(x, y)
                x, y, elevation = x.flatten(), y.flatten(), elevation.flatten()

                for dem_threshold in range(1, 4):
                    dem_pd = pd.DataFrame.from_dict({'elevation': elevation, 'x': x, 'y': y})
                    dem_pd = dem_pd[dem_pd['elevation'] == dem_threshold]
                    dem_vector = gpd.GeoDataFrame(
                        geometry=gpd.GeoSeries.from_xy(dem_pd['x'], dem_pd['y'], crs=CRS.from_string('EPSG:4326')))

                    # dem_vector = dem_vector.buffer(min(self.x_pixel_size,self.y_pixel_size),join_style=3)
                    dem_vector = dem_vector.buffer(min(self.x_pixel_size, self.y_pixel_size) / 2, cap_style=3)
                    dem_vector = dem_vector.to_crs('EPSG:4326')
                    geom_arr = []
                    # Converting GeoSeries to list of geometries
                    geoms = list(dem_vector)

                    # Converting geometries list to nested list of geometries
                    for i in range(0, len(geoms), 10000):
                        geom_arr.append(geoms[i:i + 10000])

                    # Creating multiprocessing pool to perform union operation of chunks of geometries
                    with Pool(cpu_count()) as p:
                        geom_union = p.map(union, geom_arr)

                    # Perform union operation on returned unioned geometries
                    total_union = unary_union(geom_union)

                    # Creating GeoDataFrame for total_union
                    union_vector_gdf = gpd.GeoDataFrame(geometry=gpd.GeoSeries(total_union))
                    union_vector_gdf = union_vector_gdf.simplify(0.00002)
                    union_vector_gdf = union_vector_gdf.buffer(0.00002, join_style=1).buffer(-0.00002, join_style=1)

                    union_vector_gdf = gpd.GeoDataFrame(geometry=gpd.GeoSeries(union_vector_gdf))
                    union_vector_gdf['Severity'] = dem_threshold

                    # Saving GeoDataFrame to shapefile
                    union_vector_gdf.to_file(os.path.join(self.filename, time.strftime("%Y%m%d_%H%M%S") + "_" +
                                                          VariableDefinition.VARIABLE_NAME[
                                                              'sea_water_column_thickness'] + "_severity-" + str(
                        dem_threshold) + "_rank-"+str(self.coverage.rank)+".shp"), crs='EPSG:4326')

            #     if self.coverage.rank != 0:
            #         self.coverage.comm.Send(np.ascontiguousarray(local_data), dest=0)
            #     else:
            #         # Pour le proc n°1
            #         global_data[self.coverage.map_mpi[self.coverage.rank]["dst_global_t"].start + time_index,
            #                     self.coverage.map_mpi[self.coverage.rank]["dst_global_y"],
            #                     self.coverage.map_mpi[self.coverage.rank]["dst_global_x"]] = local_data
            #
            #         # Pour les autres
            #         for source in range(1, self.coverage.size):
            #             recvbuf = np.empty([self.coverage.map_mpi[source]["dst_local_y_size"],
            #                                 self.coverage.map_mpi[source]["dst_local_x_size"]])
            #             self.coverage.comm.Recv(recvbuf, source=source)
            #
            #             global_data[self.coverage.map_mpi[source]["dst_global_t"].start + time_index,
            #                         self.coverage.map_mpi[source]["dst_global_y"],
            #                         self.coverage.map_mpi[source]["dst_global_x"]] = recvbuf
            #
            # self.coverage.comm.barrier()
            #
            # if self.coverage.rank == 0:
            #
            #     logging.info('[DefaultWriter] Writing variable \'' + str(
            #         VariableDefinition.LONG_NAME['sea_water_column_thickness']) + '\'')
            #
            #     for time_index in range(0, self.coverage.get_t_size(type="target_global")):
            #         time = self.coverage.read_axis_t(type="target_global")[time_index]
            #
            #         logging.debug('[DefaultWriter] Writing variable \'' + str(VariableDefinition.LONG_NAME[
            #                                                                       'sea_water_column_thickness']) + '\' at time \'' + str(
            #             time) + '\'')
            #
            #         global_data[time_index][global_data[time_index] == 9.96921e+36] = 0 # Remove NaN
            #         global_data[time_index][global_data[time_index] <= 0.1] = 0
            #         global_data[time_index][(global_data[time_index] > 1.2) & (global_data[time_index] != 0)] = 4
            #         global_data[time_index][(global_data[time_index] >= 0.7) & (global_data[time_index] <= 1.2)] = 3
            #         global_data[time_index][(global_data[time_index] > 0.3) & (global_data[time_index] <= 0.7)] = 2
            #         global_data[time_index][(global_data[time_index] > 0.1) & (global_data[time_index] <= 0.3)] = 1
            #
            #         file = self.driver.CreateDataSource(
            #             os.path.join(self.filename, time.strftime("%Y%m%d_%H%M%S") + "_" +
            #                          VariableDefinition.VARIABLE_NAME[
            #                              'sea_water_column_thickness'] + ".shp"))
            #
            #         # # CRS info
            #         # proj = osr.SpatialReference()
            #         # proj.SetWellKnownGeogCS("EPSG:4326")
            #         # layer = file.CreateLayer("s26", srs=proj)
            #         # fldSeverity = ogr.FieldDefn('Severity', ogr.OFTInteger)
            #         # layer.CreateField(fldSeverity)
            #         #
            #         # sq = gdal_array.OpenArray(global_data[time_index])
            #         # proj = osr.SpatialReference()
            #         # proj.SetWellKnownGeogCS("EPSG:4326")
            #         # sq.SetProjection(proj.ExportToWkt())
            #         # sq.SetGeoTransform(self.geotransform)
            #         # gdal.Polygonize(sq.GetRasterBand(1), None, layer, -1, [], callback=None)
            #
            #         x, y, elevation = self.coverage.read_axis_x(type="target_global"), self.coverage.read_axis_y(type="target_global"),  global_data[time_index]
            #         x, y = np.meshgrid(x, y)
            #         x, y, elevation = x.flatten(), y.flatten(), elevation.flatten()
            #
            #         for dem_threshold in range(1,5):
            #             dem_pd = pd.DataFrame.from_dict({'elevation': elevation, 'x': x, 'y': y})
            #             dem_pd = dem_pd[dem_pd['elevation'] == dem_threshold]
            #             dem_vector = gpd.GeoDataFrame(
            #                 geometry=gpd.GeoSeries.from_xy(dem_pd['x'], dem_pd['y'], crs= CRS.from_string('EPSG:4326')))
            #
            #             #dem_vector = dem_vector.buffer(min(self.x_pixel_size,self.y_pixel_size),join_style=3)
            #             dem_vector = dem_vector.buffer(min(self.x_pixel_size, self.y_pixel_size)/2, cap_style=3)
            #             dem_vector = dem_vector.to_crs('EPSG:4326')
            #             geom_arr = []
            #             # Converting GeoSeries to list of geometries
            #             geoms = list(dem_vector)
            #
            #             # Converting geometries list to nested list of geometries
            #             for i in range(0, len(geoms), 10000):
            #                 geom_arr.append(geoms[i:i + 10000])
            #
            #             # Creating multiprocessing pool to perform union operation of chunks of geometries
            #             with Pool(cpu_count()) as p:
            #                geom_union = p.map(union, geom_arr)
            #
            #             # Perform union operation on returned unioned geometries
            #             total_union = unary_union(geom_union)
            #
            #             # Creating GeoDataFrame for total_union
            #             union_vector_gdf = gpd.GeoDataFrame(geometry=gpd.GeoSeries(total_union))
            #
            #             #union_vector_gdf = union_vector_gdf.buffer(-min(self.x_pixel_size,self.y_pixel_size)/2,resolution=2)
            #
            #             union_vector_gdf = union_vector_gdf.simplify(0.00001)
            #             #union_vector_gdf = union_vector_gdf.buffer(min(self.x_pixel_size, self.y_pixel_size)/2,cap_style=3).buffer(-min(self.x_pixel_size, self.y_pixel_size)/2,cap_style=3)
            #             #union_vector_gdf = union_vector_gdf.simplify(min(self.x_pixel_size, self.y_pixel_size)/2)
            #             #union_vector_gdf = union_vector_gdf.buffer(0.001, join_style=1)
            #
            #             union_vector_gdf = gpd.GeoDataFrame(geometry=gpd.GeoSeries(union_vector_gdf))
            #             union_vector_gdf['Severity'] = dem_threshold
            #
            #             # Saving GeoDataFrame to shapefile
            #             union_vector_gdf.to_file(os.path.join(self.filename, time.strftime("%Y%m%d_%H%M%S") + "_" +
            #                                                    VariableDefinition.VARIABLE_NAME[
            #                                                        'sea_water_column_thickness'] +"_severity-"+str(dem_threshold)+".shp"), crs='EPSG:4326')
        else:
            raise CoverageError("DefaultWriter",
                                "The given coverage is not an instance of 'TimeCoverage' or 'TimeLevelCoverage'")

