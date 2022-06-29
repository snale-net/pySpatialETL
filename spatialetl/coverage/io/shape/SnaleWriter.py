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
from multiprocessing import Pool, cpu_count

import geopandas as gpd
import numpy as np
import pandas as pd
from rasterio.crs import CRS
from shapely.ops import unary_union

from spatialetl.coverage.TimeCoverage import TimeCoverage
from spatialetl.coverage.TimeLevelCoverage import TimeLevelCoverage
from spatialetl.coverage.io.CoverageWriter import CoverageWriter
from spatialetl.exception.CoverageError import CoverageError
from spatialetl.utils.VariableDefinition import VariableDefinition


def union(x):
    return unary_union(x)


class SnaleWriter(CoverageWriter):

    def __init__(self, cov, myFile):
        CoverageWriter.__init__(self, cov, myFile);

        if self.coverage.is_regular_grid() == False:
            raise ValueError("This writer supports only Coverage with a regular horizontal axis.")

        if os.path.isdir(self.filename) is False:
            raise ValueError("Filename has to be a directory.")

        self.rows = self.coverage.get_x_size(type="target_global")
        self.cols = self.coverage.get_y_size(type="target_global")

        xmin = np.min(self.coverage.read_axis_x(type="target_global"))
        xmax = np.max(self.coverage.read_axis_x(type="target_global"))
        ymin = np.min(self.coverage.read_axis_y(type="target_global"))
        ymax = np.max(self.coverage.read_axis_y(type="target_global"))
        self.x_pixel_size = round((xmax - xmin) / self.coverage.get_x_size(type="target_global"), 6)
        self.y_pixel_size = round((ymax - ymin) / self.coverage.get_y_size(type="target_global"), 6)

    def close(self):
        return

    # Variables

    #################
    # HYDRO
    # Sea Surface
    #################
    def write_variable_sea_water_column_thickness(self):

        if (isinstance(self.coverage, TimeCoverage) or isinstance(self.coverage, TimeLevelCoverage)):

            for time_index in range(0, self.coverage.get_t_size()):
                time = self.coverage.read_axis_t(type="target_global")[time_index]
                local_data = self.coverage.read_variable_sea_water_column_thickness_at_time(time_index)

                local_data[local_data == 9.96921e+36] = 0  # Remove NaN
                local_data[local_data <= 0.1] = 0
                local_data[(local_data > 1.2) & (local_data != 0)] = 3
                local_data[(local_data >= 0.7) & (local_data <= 1.2)] = 2
                local_data[(local_data > 0.1) & (local_data <= 0.7)] = 1
                # local_data[(local_data > 0.1) & (local_data <= 0.3)] = 1

                x, y, elevation = self.coverage.read_axis_x(type="target_local"), self.coverage.read_axis_y(
                    type="target_local"), local_data

                x, y = np.meshgrid(x, y)
                x, y, elevation = x.flatten(), y.flatten(), elevation.flatten()

                for dem_threshold in range(1, 4):
                    dem_pd = pd.DataFrame.from_dict({'elevation': elevation, 'x': x, 'y': y})
                    dem_pd = dem_pd[dem_pd['elevation'] == dem_threshold]
                    dem_vector = gpd.GeoDataFrame(
                        geometry=gpd.GeoSeries.from_xy(dem_pd['x'], dem_pd['y'], crs=CRS.from_string('EPSG:4326')))

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
                    union_vector_gdf['Datetime'] = time.strftime("%Y-%m-%d %H:%M:%S")

                    # Saving GeoDataFrame to shapefile
                    union_vector_gdf.to_file(os.path.join(self.filename, time.strftime("%Y%m%d_%H%M%S") + "_" +
                                                          VariableDefinition.VARIABLE_NAME[
                                                              'sea_water_column_thickness'] + "_severity-" + str(
                        dem_threshold) + "_rank-" + str(self.coverage.rank) + ".shp"), crs='EPSG:4326')

        else:
            raise CoverageError("DefaultWriter",
                                "The given coverage is not an instance of 'TimeCoverage' or 'TimeLevelCoverage'")
