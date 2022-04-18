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

import numpy as np
from osgeo import gdal
from osgeo import osr

from spatialetl.coverage.TimeCoverage import TimeCoverage
from spatialetl.coverage.TimeLevelCoverage import TimeLevelCoverage
from spatialetl.coverage.io.CoverageWriter import CoverageWriter
from spatialetl.exception.CoverageError import CoverageError
from spatialetl.utils.VariableDefinition import VariableDefinition
from spatialetl.utils.logger import logging


class DefaultWriter (CoverageWriter):

    def __init__(self,cov,myFile):
        CoverageWriter.__init__(self,cov,myFile);

        if self.coverage.is_regular_grid() == False:
            raise ValueError("This writer supports only Coverage with a regular horizontal axis.")

        if os.path.isdir(self.filename) is False:
            raise ValueError("Filename has to be a directory.")

        gdal.AllRegister()
        self.driver = gdal.GetDriverByName('GTiff')

        self.rows = self.coverage.get_x_size(type="target_global")
        self.cols = self.coverage.get_y_size(type="target_global")

        xmin=np.min(self.coverage.read_axis_x(type="target_global"))
        xmax=np.max(self.coverage.read_axis_x(type="target_global"))
        ymin=np.min(self.coverage.read_axis_y(type="target_global"))
        ymax=np.max(self.coverage.read_axis_y(type="target_global"))
        x_pixel_size=round((xmax-xmin)/self.coverage.get_x_size(type="target_global"),6)
        y_pixel_size=round((ymax-ymin)/self.coverage.get_y_size(type="target_global"),6)

        self.geotransform = (xmin, x_pixel_size, 0, ymin, 0, y_pixel_size)

    def close(self):
        return

    # Variables
    def write_variable_mesh_size(self):

        if self.coverage.rank == 0:
            logging.info('[DefaultWriter] Gathering data from all processors')
            global_data = np.empty(
                [self.coverage.get_y_size(type="target_global"),
                 self.coverage.get_x_size(type="target_global")])
            global_data[:] = np.nan

        local_data = self.coverage.read_variable_mesh_size()

        if self.coverage.rank != 0:
            self.coverage.comm.Send(np.ascontiguousarray(local_data), dest=0)
        else:
            # Pour le proc n°1
            global_data[self.coverage.map_mpi[self.coverage.rank]["dst_global_y"],
                        self.coverage.map_mpi[self.coverage.rank]["dst_global_x"]] = local_data

            # Pour les autres
            for source in range(1, self.coverage.size):
                recvbuf = np.empty([self.coverage.map_mpi[source]["dst_local_y_size"],
                                    self.coverage.map_mpi[source]["dst_local_x_size"]])
                self.coverage.comm.Recv(recvbuf, source=source)

                global_data[self.coverage.map_mpi[source]["dst_global_y"],
                            self.coverage.map_mpi[source]["dst_global_x"]] = recvbuf

        self.coverage.comm.barrier()

        if self.coverage.rank == 0:
            logging.info('[DefaultWriter] Writing variable \'' + str(
                VariableDefinition.LONG_NAME['mesh_size']) + '\'')

            file = self.driver.Create(
                os.path.join(self.filename, VariableDefinition.VARIABLE_NAME['mesh_size'] + ".tiff"),
                int(self.rows), int(self.cols), 1, gdal.GDT_Float64)

            # CRS info
            proj = osr.SpatialReference()
            proj.SetWellKnownGeogCS("EPSG:4326")
            file.SetProjection(proj.ExportToWkt())
            file.SetGeoTransform(self.geotransform)
            file.GetRasterBand(1).WriteArray(global_data)

    def write_variable_2D_sea_binary_mask(self):

        if self.coverage.rank == 0:
            logging.info('[DefaultWriter] Gathering data from all processors')
            global_data = np.empty(
                [self.coverage.get_y_size(type="target_global"),
                 self.coverage.get_x_size(type="target_global")])
            global_data[:] = np.nan

        local_data = self.coverage.read_variable_2d_sea_binary_mask()

        if self.coverage.rank != 0:
            self.coverage.comm.Send(np.ascontiguousarray(local_data), dest=0)
        else:
            # Pour le proc n°1
            global_data[self.coverage.map_mpi[self.coverage.rank]["dst_global_y"],
                        self.coverage.map_mpi[self.coverage.rank]["dst_global_x"]] = local_data

            # Pour les autres
            for source in range(1, self.coverage.size):
                recvbuf = np.empty([self.coverage.map_mpi[source]["dst_local_y_size"],
                                    self.coverage.map_mpi[source]["dst_local_x_size"]])
                self.coverage.comm.Recv(recvbuf, source=source)

                global_data[self.coverage.map_mpi[source]["dst_global_y"],
                            self.coverage.map_mpi[source]["dst_global_x"]] = recvbuf

        self.coverage.comm.barrier()

        if self.coverage.rank == 0:
            logging.info('[DefaultWriter] Writing variable \'' + str(
                VariableDefinition.LONG_NAME['2d_sea_binary_mask']) + '\'')

            file = self.driver.Create(
                os.path.join(self.filename, VariableDefinition.VARIABLE_NAME['2d_sea_binary_mask'] + ".tiff"),
                int(self.rows), int(self.cols), 1, gdal.GDT_Float64)

            # CRS info
            proj = osr.SpatialReference()
            proj.SetWellKnownGeogCS("EPSG:4326")
            file.SetProjection(proj.ExportToWkt())
            file.SetGeoTransform(self.geotransform)
            file.GetRasterBand(1).WriteArray(global_data)

    def write_variable_wet_binary_mask(self):

        if (isinstance(self.coverage, TimeCoverage) or isinstance(self.coverage, TimeLevelCoverage)):

            if self.coverage.rank == 0:
                logging.info('[DefaultWriter] Gathering data from all processors')
                global_data = np.empty(
                    [self.coverage.get_t_size(type="target_global"),
                     self.coverage.get_y_size(type="target_global"),
                     self.coverage.get_x_size(type="target_global")])
                global_data[:] = np.nan

            for time_index in range(0, self.coverage.get_t_size()):

                local_data = self.coverage.read_variable_wet_binary_mask_at_time(time_index)

                if self.coverage.rank != 0:
                    self.coverage.comm.Send(np.ascontiguousarray(local_data), dest=0)
                else:
                    # Pour le proc n°1
                    global_data[self.coverage.map_mpi[self.coverage.rank]["dst_global_t"].start + time_index,
                                self.coverage.map_mpi[self.coverage.rank]["dst_global_y"],
                                self.coverage.map_mpi[self.coverage.rank]["dst_global_x"]] = local_data

                    # Pour les autres
                    for source in range(1, self.coverage.size):
                        recvbuf = np.empty([self.coverage.map_mpi[source]["dst_local_y_size"],
                                            self.coverage.map_mpi[source]["dst_local_x_size"]])
                        self.coverage.comm.Recv(recvbuf, source=source)

                        global_data[self.coverage.map_mpi[source]["dst_global_t"].start + time_index,
                                    self.coverage.map_mpi[source]["dst_global_y"],
                                    self.coverage.map_mpi[source]["dst_global_x"]] = recvbuf

            self.coverage.comm.barrier()

            if self.coverage.rank == 0:

                logging.info('[DefaultWriter] Writing variable \'' + str(
                    VariableDefinition.LONG_NAME['wet_binary_mask']) + '\'')

                for time_index in range(0, self.coverage.get_t_size(type="target_global")):
                    time = self.coverage.read_axis_t(type="target_global")[time_index]

                    logging.debug('[DefaultWriter] Writing variable \'' + str(VariableDefinition.LONG_NAME[
                                                                                  'wet_binary_mask']) + '\' at time \'' + str(
                        time) + '\'')

                    file = self.driver.Create(os.path.join(self.filename, time.strftime("%Y%m%d_%H%M%S") + "_" +
                                                           VariableDefinition.VARIABLE_NAME[
                                                               'wet_binary_mask'] + ".tiff"),
                                              int(self.rows), int(self.cols), 1, gdal.GDT_Float64)

                    # CRS info
                    proj = osr.SpatialReference()
                    proj.SetWellKnownGeogCS("EPSG:4326")
                    file.SetProjection(proj.ExportToWkt())
                    file.SetGeoTransform(self.geotransform)
                    file.GetRasterBand(1).WriteArray(global_data[time_index])
        else:
            raise CoverageError("DefaultWriter","The given coverage is not an instance of 'TimeCoverage' or 'TimeLevelCoverage'")

    def write_variable_3D_sea_binary_mask(self):

        if (isinstance(self.coverage, TimeCoverage) or isinstance(self.coverage, TimeLevelCoverage)):

            if self.coverage.rank == 0:
                logging.info('[DefaultWriter] Gathering data from all processors')
                global_data = np.empty(
                    [self.coverage.get_t_size(type="target_global"),
                     self.coverage.get_y_size(type="target_global"),
                     self.coverage.get_x_size(type="target_global")])
                global_data[:] = np.nan

            for time_index in range(0, self.coverage.get_t_size()):

                local_data = self.coverage.read_variable_2D_sea_binary_mask_at_time(time_index)

                if self.coverage.rank != 0:
                    self.coverage.comm.Send(np.ascontiguousarray(local_data), dest=0)
                else:
                    # Pour le proc n°1
                    global_data[self.coverage.map_mpi[self.coverage.rank]["dst_global_t"].start + time_index,
                                self.coverage.map_mpi[self.coverage.rank]["dst_global_y"],
                                self.coverage.map_mpi[self.coverage.rank]["dst_global_x"]] = local_data

                    # Pour les autres
                    for source in range(1, self.coverage.size):
                        recvbuf = np.empty([self.coverage.map_mpi[source]["dst_local_y_size"],
                                            self.coverage.map_mpi[source]["dst_local_x_size"]])
                        self.coverage.comm.Recv(recvbuf, source=source)

                        global_data[self.coverage.map_mpi[source]["dst_global_t"].start + time_index,
                                    self.coverage.map_mpi[source]["dst_global_y"],
                                    self.coverage.map_mpi[source]["dst_global_x"]] = recvbuf

            self.coverage.comm.barrier()

            if self.coverage.rank == 0:

                logging.info('[DefaultWriter] Writing variable \'' + str(
                    VariableDefinition.LONG_NAME['2d_sea_binary_mask']) + '\'')

                for time_index in range(0, self.coverage.get_t_size(type="target_global")):
                    time = self.coverage.read_axis_t(type="target_global")[time_index]

                    logging.debug('[DefaultWriter] Writing variable \'' + str(VariableDefinition.LONG_NAME[
                                                                                  '2d_sea_binary_mask']) + '\' at time \'' + str(
                        time) + '\'')

                    file = self.driver.Create(os.path.join(self.filename, time.strftime("%Y%m%d_%H%M%S") + "_" +
                                                           VariableDefinition.VARIABLE_NAME[
                                                               '2d_sea_binary_mask'] + ".tiff"),
                                              int(self.rows), int(self.cols), 1, gdal.GDT_Float64)

                    # CRS info
                    proj = osr.SpatialReference()
                    proj.SetWellKnownGeogCS("EPSG:4326")
                    file.SetProjection(proj.ExportToWkt())
                    file.SetGeoTransform(self.geotransform)
                    file.GetRasterBand(1).WriteArray(global_data[time_index])
        else:
            raise CoverageError("DefaultWriter",
                                "The given coverage is not an instance of 'TimeCoverage' or 'TimeLevelCoverage'")

    def write_variable_3D_land_binary_mask(self):

        if (isinstance(self.coverage, TimeCoverage) or isinstance(self.coverage, TimeLevelCoverage)):

            if self.coverage.rank == 0:
                logging.info('[DefaultWriter] Gathering data from all processors')
                global_data = np.empty(
                    [self.coverage.get_t_size(type="target_global"),
                     self.coverage.get_y_size(type="target_global"),
                     self.coverage.get_x_size(type="target_global")])
                global_data[:] = np.nan

            for time_index in range(0, self.coverage.get_t_size()):

                local_data = self.coverage.read_variable_2D_land_binary_mask_at_time(time_index)

                if self.coverage.rank != 0:
                    self.coverage.comm.Send(np.ascontiguousarray(local_data), dest=0)
                else:
                    # Pour le proc n°1
                    global_data[self.coverage.map_mpi[self.coverage.rank]["dst_global_t"].start + time_index,
                                self.coverage.map_mpi[self.coverage.rank]["dst_global_y"],
                                self.coverage.map_mpi[self.coverage.rank]["dst_global_x"]] = local_data

                    # Pour les autres
                    for source in range(1, self.coverage.size):
                        recvbuf = np.empty([self.coverage.map_mpi[source]["dst_local_y_size"],
                                            self.coverage.map_mpi[source]["dst_local_x_size"]])
                        self.coverage.comm.Recv(recvbuf, source=source)

                        global_data[self.coverage.map_mpi[source]["dst_global_t"].start + time_index,
                                    self.coverage.map_mpi[source]["dst_global_y"],
                                    self.coverage.map_mpi[source]["dst_global_x"]] = recvbuf

            self.coverage.comm.barrier()

            if self.coverage.rank == 0:

                logging.info('[DefaultWriter] Writing variable \'' + str(
                    VariableDefinition.LONG_NAME['2d_land_binary_mask']) + '\'')

                for time_index in range(0, self.coverage.get_t_size(type="target_global")):
                    time = self.coverage.read_axis_t(type="target_global")[time_index]

                    logging.debug('[DefaultWriter] Writing variable \'' + str(VariableDefinition.LONG_NAME[
                                                                                  '2d_land_binary_mask']) + '\' at time \'' + str(
                        time) + '\'')

                    file = self.driver.Create(os.path.join(self.filename, time.strftime("%Y%m%d_%H%M%S") + "_" +
                                                           VariableDefinition.VARIABLE_NAME[
                                                               '2d_land_binary_mask'] + ".tiff"),
                                              int(self.rows), int(self.cols), 1, gdal.GDT_Float64)

                    # CRS info
                    proj = osr.SpatialReference()
                    proj.SetWellKnownGeogCS("EPSG:4326")
                    file.SetProjection(proj.ExportToWkt())
                    file.SetGeoTransform(self.geotransform)
                    file.GetRasterBand(1).WriteArray(global_data[time_index])
        else:
            raise CoverageError("DefaultWriter",
                                "The given coverage is not an instance of 'TimeCoverage' or 'TimeLevelCoverage'")

    #################
    # HYDRO
    # 2D
    #################
    def write_variable_bathymetry(self):

        if self.coverage.rank == 0:
            logging.info('[DefaultWriter] Gathering data from all processors')
            global_data = np.empty(
                [self.coverage.get_y_size(type="target_global"),
                 self.coverage.get_x_size(type="target_global")])
            global_data[:] = np.nan

        local_data = self.coverage.read_variable_bathymetry()

        if self.coverage.rank != 0:
            self.coverage.comm.Send(np.ascontiguousarray(local_data), dest=0)
        else:
            # Pour le proc n°1
            global_data[self.coverage.map_mpi[self.coverage.rank]["dst_global_y"],
                        self.coverage.map_mpi[self.coverage.rank]["dst_global_x"]] = local_data

            # Pour les autres
            for source in range(1, self.coverage.size):
                recvbuf = np.empty([self.coverage.map_mpi[source]["dst_local_y_size"],
                                    self.coverage.map_mpi[source]["dst_local_x_size"]])
                self.coverage.comm.Recv(recvbuf, source=source)

                global_data[self.coverage.map_mpi[source]["dst_global_y"],
                            self.coverage.map_mpi[source]["dst_global_x"]] = recvbuf

        self.coverage.comm.barrier()

        if self.coverage.rank == 0:

            logging.info('[DefaultWriter] Writing variable \'' + str(
                VariableDefinition.LONG_NAME['bathymetry']) + '\'')

            file = self.driver.Create(os.path.join(self.filename, VariableDefinition.VARIABLE_NAME['bathymetry'] + ".tiff"),
                                      int(self.rows), int(self.cols), 1, gdal.GDT_Float64)

            # CRS info
            proj = osr.SpatialReference()
            proj.SetWellKnownGeogCS("EPSG:4326")
            file.SetProjection(proj.ExportToWkt())
            file.SetGeoTransform(self.geotransform)
            file.GetRasterBand(1).WriteArray(global_data)

    def write_variable_barotropic_sea_water_velocity(self):

        if (isinstance(self.coverage, TimeCoverage) or isinstance(self.coverage, TimeLevelCoverage)):

            if self.coverage.rank == 0:
                logging.info('[DefaultWriter] Gathering data from all processors')
                global_data = np.empty(
                    [self.coverage.get_t_size(type="target_global"),
                     2,
                     self.coverage.get_y_size(type="target_global"),
                     self.coverage.get_x_size(type="target_global")])
                global_data[:] = np.nan

            for time_index in range(0, self.coverage.get_t_size()):

                local_data = self.coverage.read_variable_barotropic_sea_water_velocity_at_time(time_index)

                if self.coverage.rank != 0:
                    self.coverage.comm.Send(np.ascontiguousarray(local_data), dest=0)
                else:
                    # Pour le proc n°1
                    global_data[self.coverage.map_mpi[self.coverage.rank]["dst_global_t"].start + time_index,
                                0,
                                self.coverage.map_mpi[self.coverage.rank]["dst_global_y"],
                                self.coverage.map_mpi[self.coverage.rank]["dst_global_x"]] = local_data[0]

                    global_data[self.coverage.map_mpi[self.coverage.rank]["dst_global_t"].start + time_index,
                                1,
                                self.coverage.map_mpi[self.coverage.rank]["dst_global_y"],
                                self.coverage.map_mpi[self.coverage.rank]["dst_global_x"]] = local_data[1]

                    # Pour les autres
                    for source in range(1, self.coverage.size):
                        recvbuf = np.empty([2, self.coverage.map_mpi[source]["dst_local_y_size"],
                                            self.coverage.map_mpi[source]["dst_local_x_size"]])
                        self.coverage.comm.Recv(recvbuf, source=source)

                        global_data[self.coverage.map_mpi[source]["dst_global_t"].start + time_index,
                                    0,
                                    self.coverage.map_mpi[source]["dst_global_y"],
                                    self.coverage.map_mpi[source]["dst_global_x"]] = recvbuf[0]

                        global_data[self.coverage.map_mpi[source]["dst_global_t"].start + time_index,
                                    1,
                                    self.coverage.map_mpi[source]["dst_global_y"],
                                    self.coverage.map_mpi[source]["dst_global_x"]] = recvbuf[1]

            self.coverage.comm.barrier()

            if self.coverage.rank == 0:

                logging.info('[DefaultWriter] Writing variable \'Barotropic Sea Water Velocity\'\'')

                file = open(self.filename, "w")

                for time_index in range(0, self.coverage.get_t_size(type="target_global")):
                    time = self.coverage.read_axis_t(type="target_global")[time_index]

                    logging.info(
                        '[DefaultWriter] Writing variable \'Barotropic Sea Water Velocity\' at time \'' + str(
                            time) + '\'')
                    for vector in range(0, 2):
                        if vector == 0:

                            file = self.driver.Create(os.path.join(self.filename, time.strftime("%Y%m%d_%H%M%S") + "_" +
                                                                   VariableDefinition.VARIABLE_NAME[
                                                                       'barotropic_eastward_sea_water_velocity'] + ".tiff"),
                                                      int(self.rows), int(self.cols), 1, gdal.GDT_Float64)

                            # CRS info
                            proj = osr.SpatialReference()
                            proj.SetWellKnownGeogCS("EPSG:4326")
                            file.SetProjection(proj.ExportToWkt())
                            file.SetGeoTransform(self.geotransform)
                            file.GetRasterBand(1).WriteArray(global_data[time_index][0])

                        else:

                            file = self.driver.Create(os.path.join(self.filename, time.strftime("%Y%m%d_%H%M%S") + "_" +
                                                                   VariableDefinition.VARIABLE_NAME[
                                                                       'barotropic_northward_sea_water_velocity'] + ".tiff"),
                                                      int(self.rows), int(self.cols), 1, gdal.GDT_Float64)

                            # CRS info
                            proj = osr.SpatialReference()
                            proj.SetWellKnownGeogCS("EPSG:4326")
                            file.SetProjection(proj.ExportToWkt())
                            file.SetGeoTransform(self.geotransform)
                            file.GetRasterBand(1).WriteArray(global_data[time_index][1])

        else:
            raise CoverageError("DefaultWriter","The given coverage is not an instance of 'TimeCoverage' or 'TimeLevelCoverage'")

    #################
    # HYDRO
    # Sea Surface
    #################
    def write_variable_sea_surface_height_above_mean_sea_level(self):

        if (isinstance(self.coverage, TimeCoverage) or isinstance(self.coverage, TimeLevelCoverage)):

            if self.coverage.rank == 0:
                logging.info('[DefaultWriter] Gathering data from all processors')
                global_data = np.empty(
                    [self.coverage.get_t_size(type="target_global"),
                     self.coverage.get_y_size(type="target_global"),
                     self.coverage.get_x_size(type="target_global")])
                global_data[:] = np.nan

            for time_index in range(0, self.coverage.get_t_size()):

                local_data = self.coverage.read_variable_sea_surface_height_above_mean_sea_level_at_time(time_index)

                if self.coverage.rank != 0:
                    self.coverage.comm.Send(np.ascontiguousarray(local_data), dest=0)
                else:
                    # Pour le proc n°1
                    global_data[self.coverage.map_mpi[self.coverage.rank]["dst_global_t"].start + time_index,
                                self.coverage.map_mpi[self.coverage.rank]["dst_global_y"],
                                self.coverage.map_mpi[self.coverage.rank]["dst_global_x"]] = local_data

                    # Pour les autres
                    for source in range(1, self.coverage.size):
                        recvbuf = np.empty([self.coverage.map_mpi[source]["dst_local_y_size"],
                                            self.coverage.map_mpi[source]["dst_local_x_size"]])
                        self.coverage.comm.Recv(recvbuf, source=source)

                        global_data[self.coverage.map_mpi[source]["dst_global_t"].start + time_index,
                                    self.coverage.map_mpi[source]["dst_global_y"],
                                    self.coverage.map_mpi[source]["dst_global_x"]] = recvbuf

            self.coverage.comm.barrier()

            if self.coverage.rank == 0:

                logging.info('[DefaultWriter] Writing variable \''+str(VariableDefinition.LONG_NAME['sea_surface_height_above_mean_sea_level'])+'\'')

                for time_index in range(0, self.coverage.get_t_size(type="target_global")):
                    time = self.coverage.read_axis_t(type="target_global")[time_index]

                    logging.debug('[DefaultWriter] Writing variable \''+str(VariableDefinition.LONG_NAME['sea_surface_height_above_mean_sea_level'])+'\' at time \'' + str(time) + '\'')

                    file = self.driver.Create(os.path.join(self.filename,time.strftime("%Y%m%d_%H%M%S")+"_"+VariableDefinition.VARIABLE_NAME['sea_surface_height_above_mean_sea_level']+".tiff"), int(self.rows),int(self.cols), 1,gdal.GDT_Float64)

                    # CRS info
                    proj = osr.SpatialReference()
                    proj.SetWellKnownGeogCS("EPSG:4326")
                    file.SetProjection(proj.ExportToWkt())
                    file.SetGeoTransform(self.geotransform)
                    file.GetRasterBand(1).WriteArray(global_data[time_index])

        else:
            raise CoverageError("DefaultWriter","The given coverage is not an instance of 'TimeCoverage' or 'TimeLevelCoverage'")

    def write_variable_sea_surface_height_above_geoid(self):

        if (isinstance(self.coverage, TimeCoverage) or isinstance(self.coverage, TimeLevelCoverage)):

            if self.coverage.rank == 0:
                logging.info('[DefaultWriter] Gathering data from all processors')
                global_data = np.empty(
                    [self.coverage.get_t_size(type="target_global"),
                     self.coverage.get_y_size(type="target_global"),
                     self.coverage.get_x_size(type="target_global")])
                global_data[:] = np.nan

            for time_index in range(0, self.coverage.get_t_size()):

                local_data = self.coverage.read_variable_sea_surface_height_above_geoid_at_time(time_index)

                if self.coverage.rank != 0:
                    self.coverage.comm.Send(np.ascontiguousarray(local_data), dest=0)
                else:
                    # Pour le proc n°1
                    global_data[self.coverage.map_mpi[self.coverage.rank]["dst_global_t"].start + time_index,
                                self.coverage.map_mpi[self.coverage.rank]["dst_global_y"],
                                self.coverage.map_mpi[self.coverage.rank]["dst_global_x"]] = local_data

                    # Pour les autres
                    for source in range(1, self.coverage.size):
                        recvbuf = np.empty([self.coverage.map_mpi[source]["dst_local_y_size"],
                                            self.coverage.map_mpi[source]["dst_local_x_size"]])
                        self.coverage.comm.Recv(recvbuf, source=source)

                        global_data[self.coverage.map_mpi[source]["dst_global_t"].start + time_index,
                                    self.coverage.map_mpi[source]["dst_global_y"],
                                    self.coverage.map_mpi[source]["dst_global_x"]] = recvbuf

            self.coverage.comm.barrier()

            if self.coverage.rank == 0:

                logging.info('[DefaultWriter] Writing variable \'' + str(
                    VariableDefinition.LONG_NAME['sea_surface_height_above_geoid']) + '\'')

                for time_index in range(0, self.coverage.get_t_size(type="target_global")):
                    time = self.coverage.read_axis_t(type="target_global")[time_index]

                    logging.debug('[DefaultWriter] Writing variable \'' + str(VariableDefinition.LONG_NAME[
                                                                                  'sea_surface_height_above_geoid']) + '\' at time \'' + str(
                        time) + '\'')

                    file = self.driver.Create(os.path.join(self.filename, time.strftime("%Y%m%d_%H%M%S") + "_" +
                                                           VariableDefinition.VARIABLE_NAME[
                                                               'sea_surface_height_above_geoid'] + ".tiff"),
                                              int(self.rows), int(self.cols), 1, gdal.GDT_Float64)

                    # CRS info
                    proj = osr.SpatialReference()
                    proj.SetWellKnownGeogCS("EPSG:4326")
                    file.SetProjection(proj.ExportToWkt())
                    file.SetGeoTransform(self.geotransform)
                    file.GetRasterBand(1).WriteArray(global_data[time_index])
        else:
            raise CoverageError("DefaultWriter","The given coverage is not an instance of 'TimeCoverage' or 'TimeLevelCoverage'")

    def write_variable_sea_water_column_thickness(self):

        if (isinstance(self.coverage, TimeCoverage) or isinstance(self.coverage, TimeLevelCoverage)):

            if self.coverage.rank == 0:
                logging.info('[DefaultWriter] Gathering data from all processors')
                global_data = np.empty(
                    [self.coverage.get_t_size(type="target_global"),
                     self.coverage.get_y_size(type="target_global"),
                     self.coverage.get_x_size(type="target_global")])
                global_data[:] = np.nan

            for time_index in range(0, self.coverage.get_t_size()):

                local_data = self.coverage.read_variable_sea_water_column_thickness_at_time(time_index)

                if self.coverage.rank != 0:
                    self.coverage.comm.Send(np.ascontiguousarray(local_data), dest=0)
                else:
                    # Pour le proc n°1
                    global_data[self.coverage.map_mpi[self.coverage.rank]["dst_global_t"].start + time_index,
                                self.coverage.map_mpi[self.coverage.rank]["dst_global_y"],
                                self.coverage.map_mpi[self.coverage.rank]["dst_global_x"]] = local_data

                    # Pour les autres
                    for source in range(1, self.coverage.size):
                        recvbuf = np.empty([self.coverage.map_mpi[source]["dst_local_y_size"],
                                            self.coverage.map_mpi[source]["dst_local_x_size"]])
                        self.coverage.comm.Recv(recvbuf, source=source)

                        global_data[self.coverage.map_mpi[source]["dst_global_t"].start + time_index,
                                    self.coverage.map_mpi[source]["dst_global_y"],
                                    self.coverage.map_mpi[source]["dst_global_x"]] = recvbuf

            self.coverage.comm.barrier()

            if self.coverage.rank == 0:

                logging.info('[DefaultWriter] Writing variable \'' + str(
                    VariableDefinition.LONG_NAME['sea_water_column_thickness']) + '\'')

                for time_index in range(0, self.coverage.get_t_size(type="target_global")):
                    time = self.coverage.read_axis_t(type="target_global")[time_index]

                    logging.debug('[DefaultWriter] Writing variable \'' + str(VariableDefinition.LONG_NAME[
                                                                                  'sea_water_column_thickness']) + '\' at time \'' + str(
                        time) + '\'')

                    file = self.driver.Create(os.path.join(self.filename, time.strftime("%Y%m%d_%H%M%S") + "_" +
                                                           VariableDefinition.VARIABLE_NAME[
                                                               'sea_water_column_thickness'] + ".tiff"),
                                              int(self.rows), int(self.cols), 1, gdal.GDT_Float64)

                    # CRS info
                    proj = osr.SpatialReference()
                    proj.SetWellKnownGeogCS("EPSG:4326")
                    file.SetProjection(proj.ExportToWkt())
                    file.SetGeoTransform(self.geotransform)
                    file.GetRasterBand(1).WriteArray(global_data[time_index])
        else:
            raise CoverageError("DefaultWriter",
                                "The given coverage is not an instance of 'TimeCoverage' or 'TimeLevelCoverage'")

    def write_variable_sea_surface_temperature(self):

        if (isinstance(self.coverage, TimeCoverage) or isinstance(self.coverage, TimeLevelCoverage)):

            if self.coverage.rank == 0:
                logging.info('[DefaultWriter] Gathering data from all processors')
                global_data = np.empty(
                    [self.coverage.get_t_size(type="target_global"),
                     self.coverage.get_y_size(type="target_global"),
                     self.coverage.get_x_size(type="target_global")])
                global_data[:] = np.nan

            for time_index in range(0, self.coverage.get_t_size()):

                local_data = self.coverage.read_variable_sea_surface_temperature_at_time(time_index)

                if self.coverage.rank != 0:
                    self.coverage.comm.Send(np.ascontiguousarray(local_data), dest=0)
                else:
                    # Pour le proc n°1
                    global_data[self.coverage.map_mpi[self.coverage.rank]["dst_global_t"].start + time_index,
                                self.coverage.map_mpi[self.coverage.rank]["dst_global_y"],
                                self.coverage.map_mpi[self.coverage.rank]["dst_global_x"]] = local_data

                    # Pour les autres
                    for source in range(1, self.coverage.size):
                        recvbuf = np.empty([self.coverage.map_mpi[source]["dst_local_y_size"],
                                            self.coverage.map_mpi[source]["dst_local_x_size"]])
                        self.coverage.comm.Recv(recvbuf, source=source)

                        global_data[self.coverage.map_mpi[source]["dst_global_t"].start + time_index,
                                    self.coverage.map_mpi[source]["dst_global_y"],
                                    self.coverage.map_mpi[source]["dst_global_x"]] = recvbuf

            self.coverage.comm.barrier()

            if self.coverage.rank == 0:

                logging.info('[DefaultWriter] Writing variable \'' + str(
                    VariableDefinition.LONG_NAME['sea_surface_temperature']) + '\'')

                for time_index in range(0, self.coverage.get_t_size(type="target_global")):
                    time = self.coverage.read_axis_t(type="target_global")[time_index]

                    logging.debug('[DefaultWriter] Writing variable \'' + str(VariableDefinition.LONG_NAME[
                                                                                  'sea_surface_temperature']) + '\' at time \'' + str(
                        time) + '\'')

                    file = self.driver.Create(os.path.join(self.filename, time.strftime("%Y%m%d_%H%M%S") + "_" +
                                                           VariableDefinition.VARIABLE_NAME[
                                                               'sea_surface_temperature'] + ".tiff"),
                                              int(self.rows), int(self.cols), 1, gdal.GDT_Float64)

                    # CRS info
                    proj = osr.SpatialReference()
                    proj.SetWellKnownGeogCS("EPSG:4326")
                    file.SetProjection(proj.ExportToWkt())
                    file.SetGeoTransform(self.geotransform)
                    file.GetRasterBand(1).WriteArray(global_data[time_index])
        else:
            raise CoverageError("DefaultWriter","The given coverage is not an instance of 'TimeCoverage' or 'TimeLevelCoverage'")

    def write_variable_sea_surface_salinity(self):

        if (isinstance(self.coverage, TimeCoverage) or isinstance(self.coverage, TimeLevelCoverage)):

            if self.coverage.rank == 0:
                logging.info('[DefaultWriter] Gathering data from all processors')
                global_data = np.empty(
                    [self.coverage.get_t_size(type="target_global"),
                     self.coverage.get_y_size(type="target_global"),
                     self.coverage.get_x_size(type="target_global")])
                global_data[:] = np.nan

            for time_index in range(0, self.coverage.get_t_size()):

                local_data = self.coverage.read_variable_sea_surface_salinity_at_time(time_index)

                if self.coverage.rank != 0:
                    self.coverage.comm.Send(np.ascontiguousarray(local_data), dest=0)
                else:
                    # Pour le proc n°1
                    global_data[self.coverage.map_mpi[self.coverage.rank]["dst_global_t"].start + time_index,
                                self.coverage.map_mpi[self.coverage.rank]["dst_global_y"],
                                self.coverage.map_mpi[self.coverage.rank]["dst_global_x"]] = local_data

                    # Pour les autres
                    for source in range(1, self.coverage.size):
                        recvbuf = np.empty([self.coverage.map_mpi[source]["dst_local_y_size"],
                                            self.coverage.map_mpi[source]["dst_local_x_size"]])
                        self.coverage.comm.Recv(recvbuf, source=source)

                        global_data[self.coverage.map_mpi[source]["dst_global_t"].start + time_index,
                                    self.coverage.map_mpi[source]["dst_global_y"],
                                    self.coverage.map_mpi[source]["dst_global_x"]] = recvbuf

            self.coverage.comm.barrier()

            if self.coverage.rank == 0:

                logging.info('[DefaultWriter] Writing variable \'' + str(
                    VariableDefinition.LONG_NAME['sea_surface_salinity']) + '\'')

                for time_index in range(0, self.coverage.get_t_size(type="target_global")):
                    time = self.coverage.read_axis_t(type="target_global")[time_index]

                    logging.debug('[DefaultWriter] Writing variable \'' + str(VariableDefinition.LONG_NAME[
                                                                                  'sea_surface_salinity']) + '\' at time \'' + str(
                        time) + '\'')

                    file = self.driver.Create(os.path.join(self.filename, time.strftime("%Y%m%d_%H%M%S") + "_" +
                                                           VariableDefinition.VARIABLE_NAME[
                                                               'sea_surface_salinity'] + ".tiff"),
                                              int(self.rows), int(self.cols), 1, gdal.GDT_Float64)

                    # CRS info
                    proj = osr.SpatialReference()
                    proj.SetWellKnownGeogCS("EPSG:4326")
                    file.SetProjection(proj.ExportToWkt())
                    file.SetGeoTransform(self.geotransform)
                    file.GetRasterBand(1).WriteArray(global_data[time_index])
        else:
            raise CoverageError("DefaultWriter","The given coverage is not an instance of 'TimeCoverage' or 'TimeLevelCoverage'")

    def write_variable_sea_water_velocity_at_sea_water_surface(self):

        if (isinstance(self.coverage, TimeCoverage) or isinstance(self.coverage, TimeLevelCoverage)):

            if self.coverage.rank == 0:
                logging.info('[DefaultWriter] Gathering data from all processors')
                global_data = np.empty(
                    [self.coverage.get_t_size(type="target_global"),
                     2,
                     self.coverage.get_y_size(type="target_global"),
                     self.coverage.get_x_size(type="target_global")])
                global_data[:] = np.nan

            for time_index in range(0, self.coverage.get_t_size()):

                local_data = self.coverage.read_variable_sea_water_velocity_at_sea_water_surface_at_time(time_index)

                if self.coverage.rank != 0:
                    self.coverage.comm.Send(np.ascontiguousarray(local_data), dest=0)
                else:
                    # Pour le proc n°1
                    global_data[self.coverage.map_mpi[self.coverage.rank]["dst_global_t"].start + time_index,
                                0,
                                self.coverage.map_mpi[self.coverage.rank]["dst_global_y"],
                                self.coverage.map_mpi[self.coverage.rank]["dst_global_x"]] = local_data[0]

                    global_data[self.coverage.map_mpi[self.coverage.rank]["dst_global_t"].start + time_index,
                                1,
                                self.coverage.map_mpi[self.coverage.rank]["dst_global_y"],
                                self.coverage.map_mpi[self.coverage.rank]["dst_global_x"]] = local_data[1]

                    # Pour les autres
                    for source in range(1, self.coverage.size):
                        recvbuf = np.empty([2, self.coverage.map_mpi[source]["dst_local_y_size"],
                                            self.coverage.map_mpi[source]["dst_local_x_size"]])
                        self.coverage.comm.Recv(recvbuf, source=source)

                        global_data[self.coverage.map_mpi[source]["dst_global_t"].start + time_index,
                                    0,
                                    self.coverage.map_mpi[source]["dst_global_y"],
                                    self.coverage.map_mpi[source]["dst_global_x"]] = recvbuf[0]

                        global_data[self.coverage.map_mpi[source]["dst_global_t"].start + time_index,
                                    1,
                                    self.coverage.map_mpi[source]["dst_global_y"],
                                    self.coverage.map_mpi[source]["dst_global_x"]] = recvbuf[1]

            self.coverage.comm.barrier()

            if self.coverage.rank == 0:

                logging.info('[DefaultWriter] Writing variable \'Sea Water Velocity at Sea Water Surface\'\'')

                file = open(self.filename, "w")

                for time_index in range(0, self.coverage.get_t_size(type="target_global")):
                    time = self.coverage.read_axis_t(type="target_global")[time_index]

                    logging.info(
                        '[DefaultWriter] Writing variable \'Sea Water Velocity at Sea Water Surface\' at time \'' + str(
                            time) + '\'')
                    for vector in range(0, 2):
                        if vector == 0:

                            file = self.driver.Create(os.path.join(self.filename, time.strftime("%Y%m%d_%H%M%S") + "_" +
                                                                   VariableDefinition.VARIABLE_NAME[
                                                                       'eastward_sea_water_velocity_at_sea_water_surface'] + ".tiff"),
                                                      int(self.rows), int(self.cols), 1, gdal.GDT_Float64)

                            # CRS info
                            proj = osr.SpatialReference()
                            proj.SetWellKnownGeogCS("EPSG:4326")
                            file.SetProjection(proj.ExportToWkt())
                            file.SetGeoTransform(self.geotransform)
                            file.GetRasterBand(1).WriteArray(global_data[time_index][0])

                        else:

                            file = self.driver.Create(os.path.join(self.filename, time.strftime("%Y%m%d_%H%M%S") + "_" +
                                                                   VariableDefinition.VARIABLE_NAME[
                                                                       'northward_sea_water_velocity_at_sea_water_surface'] + ".tiff"),
                                                      int(self.rows), int(self.cols), 1, gdal.GDT_Float64)

                            # CRS info
                            proj = osr.SpatialReference()
                            proj.SetWellKnownGeogCS("EPSG:4326")
                            file.SetProjection(proj.ExportToWkt())
                            file.SetGeoTransform(self.geotransform)
                            file.GetRasterBand(1).WriteArray(global_data[time_index][1])

        else:
            raise CoverageError("DefaultWriter","The given coverage is not an instance of 'TimeCoverage' or 'TimeLevelCoverage'")

    #################
    # HYDRO
    # Ground level
    #################

    def write_variable_sea_water_temperature_at_ground_level(self):

        if (isinstance(self.coverage, TimeCoverage) or isinstance(self.coverage, TimeLevelCoverage)):

            if self.coverage.rank == 0:
                logging.info('[DefaultWriter] Gathering data from all processors')
                global_data = np.empty(
                    [self.coverage.get_t_size(type="target_global"),
                     self.coverage.get_y_size(type="target_global"),
                     self.coverage.get_x_size(type="target_global")])
                global_data[:] = np.nan

            for time_index in range(0, self.coverage.get_t_size()):

                local_data = self.coverage.read_variable_sea_water_temperature_at_ground_level_at_time(time_index)

                if self.coverage.rank != 0:
                    self.coverage.comm.Send(np.ascontiguousarray(local_data), dest=0)
                else:
                    # Pour le proc n°1
                    global_data[self.coverage.map_mpi[self.coverage.rank]["dst_global_t"].start + time_index,
                                self.coverage.map_mpi[self.coverage.rank]["dst_global_y"],
                                self.coverage.map_mpi[self.coverage.rank]["dst_global_x"]] = local_data

                    # Pour les autres
                    for source in range(1, self.coverage.size):
                        recvbuf = np.empty([self.coverage.map_mpi[source]["dst_local_y_size"],
                                            self.coverage.map_mpi[source]["dst_local_x_size"]])
                        self.coverage.comm.Recv(recvbuf, source=source)

                        global_data[self.coverage.map_mpi[source]["dst_global_t"].start + time_index,
                                    self.coverage.map_mpi[source]["dst_global_y"],
                                    self.coverage.map_mpi[source]["dst_global_x"]] = recvbuf

            self.coverage.comm.barrier()

            if self.coverage.rank == 0:

                logging.info('[DefaultWriter] Writing variable \'' + str(
                    VariableDefinition.LONG_NAME['sea_water_temperature_at_ground_level']) + '\'')

                for time_index in range(0, self.coverage.get_t_size(type="target_global")):
                    time = self.coverage.read_axis_t(type="target_global")[time_index]

                    logging.debug('[DefaultWriter] Writing variable \'' + str(VariableDefinition.LONG_NAME[
                                                                                  'sea_water_temperature_at_ground_level']) + '\' at time \'' + str(
                        time) + '\'')

                    file = self.driver.Create(os.path.join(self.filename, time.strftime("%Y%m%d_%H%M%S") + "_" +
                                                           VariableDefinition.VARIABLE_NAME[
                                                               'sea_water_temperature_at_ground_level'] + ".tiff"),
                                              int(self.rows), int(self.cols), 1, gdal.GDT_Float64)

                    # CRS info
                    proj = osr.SpatialReference()
                    proj.SetWellKnownGeogCS("EPSG:4326")
                    file.SetProjection(proj.ExportToWkt())
                    file.SetGeoTransform(self.geotransform)
                    file.GetRasterBand(1).WriteArray(global_data[time_index])
        else:
            raise CoverageError("DefaultWriter","The given coverage is not an instance of 'TimeCoverage' or 'TimeLevelCoverage'")

    def write_variable_sea_water_salinity_at_ground_level(self):

        if (isinstance(self.coverage, TimeCoverage) or isinstance(self.coverage, TimeLevelCoverage)):

            if self.coverage.rank == 0:
                logging.info('[DefaultWriter] Gathering data from all processors')
                global_data = np.empty(
                    [self.coverage.get_t_size(type="target_global"),
                     self.coverage.get_y_size(type="target_global"),
                     self.coverage.get_x_size(type="target_global")])
                global_data[:] = np.nan

            for time_index in range(0, self.coverage.get_t_size()):

                local_data = self.coverage.read_variable_sea_water_salinity_at_ground_level_at_time(time_index)

                if self.coverage.rank != 0:
                    self.coverage.comm.Send(np.ascontiguousarray(local_data), dest=0)
                else:
                    # Pour le proc n°1
                    global_data[self.coverage.map_mpi[self.coverage.rank]["dst_global_t"].start + time_index,
                                self.coverage.map_mpi[self.coverage.rank]["dst_global_y"],
                                self.coverage.map_mpi[self.coverage.rank]["dst_global_x"]] = local_data

                    # Pour les autres
                    for source in range(1, self.coverage.size):
                        recvbuf = np.empty([self.coverage.map_mpi[source]["dst_local_y_size"],
                                            self.coverage.map_mpi[source]["dst_local_x_size"]])
                        self.coverage.comm.Recv(recvbuf, source=source)

                        global_data[self.coverage.map_mpi[source]["dst_global_t"].start + time_index,
                                    self.coverage.map_mpi[source]["dst_global_y"],
                                    self.coverage.map_mpi[source]["dst_global_x"]] = recvbuf

            self.coverage.comm.barrier()

            if self.coverage.rank == 0:

                logging.info('[DefaultWriter] Writing variable \'' + str(
                    VariableDefinition.LONG_NAME['sea_water_salinity_at_ground_level']) + '\'')

                for time_index in range(0, self.coverage.get_t_size(type="target_global")):
                    time = self.coverage.read_axis_t(type="target_global")[time_index]

                    logging.debug('[DefaultWriter] Writing variable \'' + str(VariableDefinition.LONG_NAME[
                                                                                  'sea_water_salinity_at_ground_level']) + '\' at time \'' + str(
                        time) + '\'')

                    file = self.driver.Create(os.path.join(self.filename, time.strftime("%Y%m%d_%H%M%S") + "_" +
                                                           VariableDefinition.VARIABLE_NAME[
                                                               'sea_water_salinity_at_ground_level'] + ".tiff"),
                                              int(self.rows), int(self.cols), 1, gdal.GDT_Float64)

                    # CRS info
                    proj = osr.SpatialReference()
                    proj.SetWellKnownGeogCS("EPSG:4326")
                    file.SetProjection(proj.ExportToWkt())
                    file.SetGeoTransform(self.geotransform)
                    file.GetRasterBand(1).WriteArray(global_data[time_index])
        else:
            raise CoverageError("DefaultWriter","The given coverage is not an instance of 'TimeCoverage' or 'TimeLevelCoverage'")

    def write_variable_sea_water_velocity_at_ground_level(self):

        if (isinstance(self.coverage, TimeCoverage) or isinstance(self.coverage, TimeLevelCoverage)):

            if self.coverage.rank == 0:
                logging.info('[DefaultWriter] Gathering data from all processors')
                global_data = np.empty(
                    [self.coverage.get_t_size(type="target_global"),
                     2,
                     self.coverage.get_y_size(type="target_global"),
                     self.coverage.get_x_size(type="target_global")])
                global_data[:] = np.nan

            for time_index in range(0, self.coverage.get_t_size()):

                local_data = self.coverage.read_variable_sea_water_velocity_at_ground_level_at_time(time_index)

                if self.coverage.rank != 0:
                    self.coverage.comm.Send(np.ascontiguousarray(local_data), dest=0)
                else:
                    # Pour le proc n°1
                    global_data[self.coverage.map_mpi[self.coverage.rank]["dst_global_t"].start + time_index,
                                0,
                                self.coverage.map_mpi[self.coverage.rank]["dst_global_y"],
                                self.coverage.map_mpi[self.coverage.rank]["dst_global_x"]] = local_data[0]

                    global_data[self.coverage.map_mpi[self.coverage.rank]["dst_global_t"].start + time_index,
                                1,
                                self.coverage.map_mpi[self.coverage.rank]["dst_global_y"],
                                self.coverage.map_mpi[self.coverage.rank]["dst_global_x"]] = local_data[1]

                    # Pour les autres
                    for source in range(1, self.coverage.size):
                        recvbuf = np.empty([2, self.coverage.map_mpi[source]["dst_local_y_size"],
                                            self.coverage.map_mpi[source]["dst_local_x_size"]])
                        self.coverage.comm.Recv(recvbuf, source=source)

                        global_data[self.coverage.map_mpi[source]["dst_global_t"].start + time_index,
                                    0,
                                    self.coverage.map_mpi[source]["dst_global_y"],
                                    self.coverage.map_mpi[source]["dst_global_x"]] = recvbuf[0]

                        global_data[self.coverage.map_mpi[source]["dst_global_t"].start + time_index,
                                    1,
                                    self.coverage.map_mpi[source]["dst_global_y"],
                                    self.coverage.map_mpi[source]["dst_global_x"]] = recvbuf[1]

            self.coverage.comm.barrier()

            if self.coverage.rank == 0:

                logging.info('[DefaultWriter] Writing variable \'Sea Water Velocity at Ground Level\'\'')

                file = open(self.filename, "w")

                for time_index in range(0, self.coverage.get_t_size(type="target_global")):
                    time = self.coverage.read_axis_t(type="target_global")[time_index]

                    logging.info(
                        '[DefaultWriter] Writing variable \'Sea Water Velocity at Ground Level\' at time \'' + str(
                            time) + '\'')
                    for vector in range(0, 2):
                        if vector == 0:

                            file = self.driver.Create(os.path.join(self.filename, time.strftime("%Y%m%d_%H%M%S") + "_" +
                                                                   VariableDefinition.VARIABLE_NAME[
                                                                       'eastward_sea_water_velocity_at_ground_level'] + ".tiff"),
                                                      int(self.rows), int(self.cols), 1, gdal.GDT_Float64)

                            # CRS info
                            proj = osr.SpatialReference()
                            proj.SetWellKnownGeogCS("EPSG:4326")
                            file.SetProjection(proj.ExportToWkt())
                            file.SetGeoTransform(self.geotransform)
                            file.GetRasterBand(1).WriteArray(global_data[time_index][0])

                        else:

                            file = self.driver.Create(os.path.join(self.filename, time.strftime("%Y%m%d_%H%M%S") + "_" +
                                                                   VariableDefinition.VARIABLE_NAME[
                                                                       'northward_sea_water_velocity_at_ground_level'] + ".tiff"),
                                                      int(self.rows), int(self.cols), 1, gdal.GDT_Float64)

                            # CRS info
                            proj = osr.SpatialReference()
                            proj.SetWellKnownGeogCS("EPSG:4326")
                            file.SetProjection(proj.ExportToWkt())
                            file.SetGeoTransform(self.geotransform)
                            file.GetRasterBand(1).WriteArray(global_data[time_index][1])
        else:
            raise CoverageError("DefaultWriter","The given coverage is not an instance of 'TimeCoverage' or 'TimeLevelCoverage'")


    #################
    # HYDRO
    # 3D
    #################

    # def write_variable_depth(self):
    #
    #     if (isinstance(self.coverage, LevelCoverage) or isinstance(self.coverage, TimeLevelCoverage)):
    #
    #         #if VariableDefinition.VARIABLE_NAME['depth_sigma'] in self.ncfile.variables:
    #         #self.coverage.read_variable_depth_at_depth(level)
    #     else:
    #         raise CoverageError("DefaultWriter",
    #                             "The given coverage is not an instance of 'LevelCoverage' or 'TimeLevelCoverage'")
    #
    # def write_variable_sea_water_temperature(self):
    #
    #     if (isinstance(self.coverage, TimeLevelCoverage)):
    #
    #         #if VariableDefinition.VARIABLE_NAME['sea_water_temperature'] in self.ncfile.variables:
    #         level_index = 0
    #         for level in self.coverage.read_axis_z():
    #             # Pas d'interpolation temporelle donc on parcours les index du temps
    #             var[
    #             self.coverage.map_mpi[self.coverage.rank]["dst_global_t"].start + time_index:self.coverage.map_mpi[self.coverage.rank]["dst_global_t"].start + time_index + 1,
    #             level_index:level_index + 1,
    #             self.coverage.map_mpi[self.coverage.rank]["dst_global_y"],
    #             self.coverage.map_mpi[self.coverage.rank]["dst_global_x"]
    #             ] = self.coverage.read_variable_sea_water_temperature_at_time_and_depth(time, level)
    #
    #             level_index += 1
    #
    #     else:
    #         raise CoverageError("DefaultWriter","The given coverage is not an instance of 'TimeLevelCoverage'")
    #
    # def write_variable_sea_water_salinity(self):
    #
    #     if (isinstance(self.coverage, TimeLevelCoverage)):
    #
    #         #if VariableDefinition.VARIABLE_NAME['sea_water_salinity'] in self.ncfile.variables:
    #         level_index = 0
    #         for level in self.coverage.read_axis_z():
    #             # Pas d'interpolation temporelle donc on parcours les index du temps
    #             var[
    #             self.coverage.map_mpi[self.coverage.rank]["dst_global_t"].start + time_index:self.coverage.map_mpi[self.coverage.rank]["dst_global_t"].start + time_index + 1,
    #             level_index:level_index + 1,
    #             self.coverage.map_mpi[self.coverage.rank]["dst_global_y"],
    #             self.coverage.map_mpi[self.coverage.rank]["dst_global_x"]
    #             ] = self.coverage.read_variable_sea_water_salinity_at_time_and_depth(time, level)
    #
    #             level_index += 1
    #
    #
    #     else:
    #         raise CoverageError("DefaultWriter","The given coverage is not an instance of 'TimeLevelCoverage'")
    #
    # def write_variable_baroclinic_sea_water_velocity(self):
    #
    #     if (isinstance(self.coverage, TimeLevelCoverage)):
    #
    #         #if VariableDefinition.VARIABLE_NAME['baroclinic_eastward_sea_water_velocity'] in self.ncfile.variables:
    #         #if VariableDefinition.VARIABLE_NAME['baroclinic_northward_sea_water_velocity'] in self.ncfile.variables:
    #
    #     else:
    #         raise CoverageError("DefaultWriter","The given coverage is not an instance of 'TimeLevelCoverage'")

    #################
    # WAVES
    # Sea Surface
    #################
    def write_variable_sea_surface_wave_significant_height(self):

        if (isinstance(self.coverage, TimeCoverage) or isinstance(self.coverage, TimeLevelCoverage)):

            if self.coverage.rank == 0:
                logging.info('[DefaultWriter] Gathering data from all processors')
                global_data = np.empty(
                    [self.coverage.get_t_size(type="target_global"),
                     self.coverage.get_y_size(type="target_global"),
                     self.coverage.get_x_size(type="target_global")])
                global_data[:] = np.nan

            for time_index in range(0, self.coverage.get_t_size()):

                local_data = self.coverage.read_variable_sea_surface_wave_significant_height_at_time(time_index)

                if self.coverage.rank != 0:
                    self.coverage.comm.Send(np.ascontiguousarray(local_data), dest=0)
                else:
                    # Pour le proc n°1
                    global_data[self.coverage.map_mpi[self.coverage.rank]["dst_global_t"].start + time_index,
                                self.coverage.map_mpi[self.coverage.rank]["dst_global_y"],
                                self.coverage.map_mpi[self.coverage.rank]["dst_global_x"]] = local_data

                    # Pour les autres
                    for source in range(1, self.coverage.size):
                        recvbuf = np.empty([self.coverage.map_mpi[source]["dst_local_y_size"],
                                            self.coverage.map_mpi[source]["dst_local_x_size"]])
                        self.coverage.comm.Recv(recvbuf, source=source)

                        global_data[self.coverage.map_mpi[source]["dst_global_t"].start + time_index,
                                    self.coverage.map_mpi[source]["dst_global_y"],
                                    self.coverage.map_mpi[source]["dst_global_x"]] = recvbuf

            self.coverage.comm.barrier()

            if self.coverage.rank == 0:

                logging.info('[DefaultWriter] Writing variable \'' + str(
                    VariableDefinition.LONG_NAME['sea_surface_wave_significant_height']) + '\'')

                for time_index in range(0, self.coverage.get_t_size(type="target_global")):
                    time = self.coverage.read_axis_t(type="target_global")[time_index]

                    logging.debug('[DefaultWriter] Writing variable \'' + str(VariableDefinition.LONG_NAME[
                                                                                  'sea_surface_wave_significant_height']) + '\' at time \'' + str(
                        time) + '\'')

                    file = self.driver.Create(os.path.join(self.filename, time.strftime("%Y%m%d_%H%M%S") + "_" +
                                                           VariableDefinition.VARIABLE_NAME[
                                                               'sea_surface_wave_significant_height'] + ".tiff"),
                                              int(self.rows), int(self.cols), 1, gdal.GDT_Float64)

                    # CRS info
                    proj = osr.SpatialReference()
                    proj.SetWellKnownGeogCS("EPSG:4326")
                    file.SetProjection(proj.ExportToWkt())
                    file.SetGeoTransform(self.geotransform)
                    file.GetRasterBand(1).WriteArray(global_data[time_index])
        else:
            raise CoverageError("DefaultWriter","The given coverage is not an instance of 'TimeCoverage' or 'TimeLevelCoverage'")

    def write_variable_sea_surface_wave_breaking_height(self):

        if (isinstance(self.coverage, TimeCoverage) or isinstance(self.coverage, TimeLevelCoverage)):

            if self.coverage.rank == 0:
                logging.info('[DefaultWriter] Gathering data from all processors')
                global_data = np.empty(
                    [self.coverage.get_t_size(type="target_global"),
                     self.coverage.get_y_size(type="target_global"),
                     self.coverage.get_x_size(type="target_global")])
                global_data[:] = np.nan

            for time_index in range(0, self.coverage.get_t_size()):

                local_data = self.coverage.read_variable_sea_surface_wave_breaking_height_at_time(time_index)

                if self.coverage.rank != 0:
                    self.coverage.comm.Send(np.ascontiguousarray(local_data), dest=0)
                else:
                    # Pour le proc n°1
                    global_data[self.coverage.map_mpi[self.coverage.rank]["dst_global_t"].start + time_index,
                                self.coverage.map_mpi[self.coverage.rank]["dst_global_y"],
                                self.coverage.map_mpi[self.coverage.rank]["dst_global_x"]] = local_data

                    # Pour les autres
                    for source in range(1, self.coverage.size):
                        recvbuf = np.empty([self.coverage.map_mpi[source]["dst_local_y_size"],
                                            self.coverage.map_mpi[source]["dst_local_x_size"]])
                        self.coverage.comm.Recv(recvbuf, source=source)

                        global_data[self.coverage.map_mpi[source]["dst_global_t"].start + time_index,
                                    self.coverage.map_mpi[source]["dst_global_y"],
                                    self.coverage.map_mpi[source]["dst_global_x"]] = recvbuf

            self.coverage.comm.barrier()

            if self.coverage.rank == 0:

                logging.info('[DefaultWriter] Writing variable \'' + str(
                    VariableDefinition.LONG_NAME['sea_surface_wave_breaking_height']) + '\'')

                for time_index in range(0, self.coverage.get_t_size(type="target_global")):
                    time = self.coverage.read_axis_t(type="target_global")[time_index]

                    logging.debug('[DefaultWriter] Writing variable \'' + str(VariableDefinition.LONG_NAME[
                                                                                  'sea_surface_wave_breaking_height']) + '\' at time \'' + str(
                        time) + '\'')

                    file = self.driver.Create(os.path.join(self.filename, time.strftime("%Y%m%d_%H%M%S") + "_" +
                                                           VariableDefinition.VARIABLE_NAME[
                                                               'sea_surface_wave_breaking_height'] + ".tiff"),
                                              int(self.rows), int(self.cols), 1, gdal.GDT_Float64)

                    # CRS info
                    proj = osr.SpatialReference()
                    proj.SetWellKnownGeogCS("EPSG:4326")
                    file.SetProjection(proj.ExportToWkt())
                    file.SetGeoTransform(self.geotransform)
                    file.GetRasterBand(1).WriteArray(global_data[time_index])
        else:
            raise CoverageError("DefaultWriter","The given coverage is not an instance of 'TimeCoverage' or 'TimeLevelCoverage'")

    def write_variable_sea_surface_wave_mean_period(self):

        if (isinstance(self.coverage, TimeCoverage) or isinstance(self.coverage, TimeLevelCoverage)):

            if self.coverage.rank == 0:
                logging.info('[DefaultWriter] Gathering data from all processors')
                global_data = np.empty(
                    [self.coverage.get_t_size(type="target_global"),
                     self.coverage.get_y_size(type="target_global"),
                     self.coverage.get_x_size(type="target_global")])
                global_data[:] = np.nan

            for time_index in range(0, self.coverage.get_t_size()):

                local_data = self.coverage.read_variable_sea_surface_wave_mean_period_at_time(time_index)

                if self.coverage.rank != 0:
                    self.coverage.comm.Send(np.ascontiguousarray(local_data), dest=0)
                else:
                    # Pour le proc n°1
                    global_data[self.coverage.map_mpi[self.coverage.rank]["dst_global_t"].start + time_index,
                                self.coverage.map_mpi[self.coverage.rank]["dst_global_y"],
                                self.coverage.map_mpi[self.coverage.rank]["dst_global_x"]] = local_data

                    # Pour les autres
                    for source in range(1, self.coverage.size):
                        recvbuf = np.empty([self.coverage.map_mpi[source]["dst_local_y_size"],
                                            self.coverage.map_mpi[source]["dst_local_x_size"]])
                        self.coverage.comm.Recv(recvbuf, source=source)

                        global_data[self.coverage.map_mpi[source]["dst_global_t"].start + time_index,
                                    self.coverage.map_mpi[source]["dst_global_y"],
                                    self.coverage.map_mpi[source]["dst_global_x"]] = recvbuf

            self.coverage.comm.barrier()

            if self.coverage.rank == 0:

                logging.info('[DefaultWriter] Writing variable \'' + str(
                    VariableDefinition.LONG_NAME['sea_surface_wave_mean_period']) + '\'')

                for time_index in range(0, self.coverage.get_t_size(type="target_global")):
                    time = self.coverage.read_axis_t(type="target_global")[time_index]

                    logging.debug('[DefaultWriter] Writing variable \'' + str(VariableDefinition.LONG_NAME[
                                                                                  'sea_surface_wave_mean_period']) + '\' at time \'' + str(
                        time) + '\'')

                    file = self.driver.Create(os.path.join(self.filename, time.strftime("%Y%m%d_%H%M%S") + "_" +
                                                           VariableDefinition.VARIABLE_NAME[
                                                               'sea_surface_wave_mean_period'] + ".tiff"),
                                              int(self.rows), int(self.cols), 1, gdal.GDT_Float64)

                    # CRS info
                    proj = osr.SpatialReference()
                    proj.SetWellKnownGeogCS("EPSG:4326")
                    file.SetProjection(proj.ExportToWkt())
                    file.SetGeoTransform(self.geotransform)
                    file.GetRasterBand(1).WriteArray(global_data[time_index])
        else:
            raise CoverageError("DefaultWriter","The given coverage is not an instance of 'TimeCoverage' or 'TimeLevelCoverage'")

    def write_variable_sea_surface_wave_peak_period(self):

        if (isinstance(self.coverage, TimeCoverage) or isinstance(self.coverage, TimeLevelCoverage)):

            if self.coverage.rank == 0:
                logging.info('[DefaultWriter] Gathering data from all processors')
                global_data = np.empty(
                    [self.coverage.get_t_size(type="target_global"),
                     self.coverage.get_y_size(type="target_global"),
                     self.coverage.get_x_size(type="target_global")])
                global_data[:] = np.nan

            for time_index in range(0, self.coverage.get_t_size()):

                local_data = self.coverage.read_variable_sea_surface_wave_peak_period_at_time(time_index)

                if self.coverage.rank != 0:
                    self.coverage.comm.Send(np.ascontiguousarray(local_data), dest=0)
                else:
                    # Pour le proc n°1
                    global_data[self.coverage.map_mpi[self.coverage.rank]["dst_global_t"].start + time_index,
                                self.coverage.map_mpi[self.coverage.rank]["dst_global_y"],
                                self.coverage.map_mpi[self.coverage.rank]["dst_global_x"]] = local_data

                    # Pour les autres
                    for source in range(1, self.coverage.size):
                        recvbuf = np.empty([self.coverage.map_mpi[source]["dst_local_y_size"],
                                            self.coverage.map_mpi[source]["dst_local_x_size"]])
                        self.coverage.comm.Recv(recvbuf, source=source)

                        global_data[self.coverage.map_mpi[source]["dst_global_t"].start + time_index,
                                    self.coverage.map_mpi[source]["dst_global_y"],
                                    self.coverage.map_mpi[source]["dst_global_x"]] = recvbuf

            self.coverage.comm.barrier()

            if self.coverage.rank == 0:

                logging.info('[DefaultWriter] Writing variable \'' + str(
                    VariableDefinition.LONG_NAME['sea_surface_wave_peak_period']) + '\'')

                for time_index in range(0, self.coverage.get_t_size(type="target_global")):
                    time = self.coverage.read_axis_t(type="target_global")[time_index]

                    logging.debug('[DefaultWriter] Writing variable \'' + str(VariableDefinition.LONG_NAME[
                                                                                  'sea_surface_wave_peak_period']) + '\' at time \'' + str(
                        time) + '\'')

                    file = self.driver.Create(os.path.join(self.filename, time.strftime("%Y%m%d_%H%M%S") + "_" +
                                                           VariableDefinition.VARIABLE_NAME[
                                                               'sea_surface_wave_peak_period'] + ".tiff"),
                                              int(self.rows), int(self.cols), 1, gdal.GDT_Float64)

                    # CRS info
                    proj = osr.SpatialReference()
                    proj.SetWellKnownGeogCS("EPSG:4326")
                    file.SetProjection(proj.ExportToWkt())
                    file.SetGeoTransform(self.geotransform)
                    file.GetRasterBand(1).WriteArray(global_data[time_index])
        else:
            raise CoverageError("DefaultWriter","The given coverage is not an instance of 'TimeCoverage' or 'TimeLevelCoverage'")

    def write_variable_sea_surface_wave_from_direction(self):

        if (isinstance(self.coverage, TimeCoverage) or isinstance(self.coverage, TimeLevelCoverage)):

            if self.coverage.rank == 0:
                logging.info('[DefaultWriter] Gathering data from all processors')
                global_data = np.empty(
                    [self.coverage.get_t_size(type="target_global"),
                     self.coverage.get_y_size(type="target_global"),
                     self.coverage.get_x_size(type="target_global")])
                global_data[:] = np.nan

            for time_index in range(0, self.coverage.get_t_size()):

                local_data = self.coverage.read_variable_sea_surface_wave_from_direction_at_time(time_index)

                if self.coverage.rank != 0:
                    self.coverage.comm.Send(np.ascontiguousarray(local_data), dest=0)
                else:
                    # Pour le proc n°1
                    global_data[self.coverage.map_mpi[self.coverage.rank]["dst_global_t"].start + time_index,
                                self.coverage.map_mpi[self.coverage.rank]["dst_global_y"],
                                self.coverage.map_mpi[self.coverage.rank]["dst_global_x"]] = local_data

                    # Pour les autres
                    for source in range(1, self.coverage.size):
                        recvbuf = np.empty([self.coverage.map_mpi[source]["dst_local_y_size"],
                                            self.coverage.map_mpi[source]["dst_local_x_size"]])
                        self.coverage.comm.Recv(recvbuf, source=source)

                        global_data[self.coverage.map_mpi[source]["dst_global_t"].start + time_index,
                                    self.coverage.map_mpi[source]["dst_global_y"],
                                    self.coverage.map_mpi[source]["dst_global_x"]] = recvbuf

            self.coverage.comm.barrier()

            if self.coverage.rank == 0:

                logging.info('[DefaultWriter] Writing variable \'' + str(
                    VariableDefinition.LONG_NAME['sea_surface_wave_from_direction']) + '\'')

                for time_index in range(0, self.coverage.get_t_size(type="target_global")):
                    time = self.coverage.read_axis_t(type="target_global")[time_index]

                    logging.debug('[DefaultWriter] Writing variable \'' + str(VariableDefinition.LONG_NAME[
                                                                                  'sea_surface_wave_from_direction']) + '\' at time \'' + str(
                        time) + '\'')

                    file = self.driver.Create(os.path.join(self.filename, time.strftime("%Y%m%d_%H%M%S") + "_" +
                                                           VariableDefinition.VARIABLE_NAME[
                                                               'sea_surface_wave_from_direction'] + ".tiff"),
                                              int(self.rows), int(self.cols), 1, gdal.GDT_Float64)

                    # CRS info
                    proj = osr.SpatialReference()
                    proj.SetWellKnownGeogCS("EPSG:4326")
                    file.SetProjection(proj.ExportToWkt())
                    file.SetGeoTransform(self.geotransform)
                    file.GetRasterBand(1).WriteArray(global_data[time_index])
        else:
            raise CoverageError("DefaultWriter","The given coverage is not an instance of 'TimeCoverage' or 'TimeLevelCoverage'")

    def write_variable_sea_surface_wave_to_direction(self):

        if (isinstance(self.coverage, TimeCoverage) or isinstance(self.coverage, TimeLevelCoverage)):

            if self.coverage.rank == 0:
                logging.info('[DefaultWriter] Gathering data from all processors')
                global_data = np.empty(
                    [self.coverage.get_t_size(type="target_global"),
                     self.coverage.get_y_size(type="target_global"),
                     self.coverage.get_x_size(type="target_global")])
                global_data[:] = np.nan

            for time_index in range(0, self.coverage.get_t_size()):

                local_data = self.coverage.read_variable_sea_surface_wave_to_direction_at_time(time_index)

                if self.coverage.rank != 0:
                    self.coverage.comm.Send(np.ascontiguousarray(local_data), dest=0)
                else:
                    # Pour le proc n°1
                    global_data[self.coverage.map_mpi[self.coverage.rank]["dst_global_t"].start + time_index,
                                self.coverage.map_mpi[self.coverage.rank]["dst_global_y"],
                                self.coverage.map_mpi[self.coverage.rank]["dst_global_x"]] = local_data

                    # Pour les autres
                    for source in range(1, self.coverage.size):
                        recvbuf = np.empty([self.coverage.map_mpi[source]["dst_local_y_size"],
                                            self.coverage.map_mpi[source]["dst_local_x_size"]])
                        self.coverage.comm.Recv(recvbuf, source=source)

                        global_data[self.coverage.map_mpi[source]["dst_global_t"].start + time_index,
                                    self.coverage.map_mpi[source]["dst_global_y"],
                                    self.coverage.map_mpi[source]["dst_global_x"]] = recvbuf

            self.coverage.comm.barrier()

            if self.coverage.rank == 0:

                logging.info('[DefaultWriter] Writing variable \'' + str(
                    VariableDefinition.LONG_NAME['sea_surface_wave_to_direction']) + '\'')

                for time_index in range(0, self.coverage.get_t_size(type="target_global")):
                    time = self.coverage.read_axis_t(type="target_global")[time_index]

                    logging.debug('[DefaultWriter] Writing variable \'' + str(VariableDefinition.LONG_NAME[
                                                                                  'sea_surface_wave_to_direction']) + '\' at time \'' + str(
                        time) + '\'')

                    file = self.driver.Create(os.path.join(self.filename, time.strftime("%Y%m%d_%H%M%S") + "_" +
                                                           VariableDefinition.VARIABLE_NAME[
                                                               'sea_surface_wave_to_direction'] + ".tiff"),
                                              int(self.rows), int(self.cols), 1, gdal.GDT_Float64)

                    # CRS info
                    proj = osr.SpatialReference()
                    proj.SetWellKnownGeogCS("EPSG:4326")
                    file.SetProjection(proj.ExportToWkt())
                    file.SetGeoTransform(self.geotransform)
                    file.GetRasterBand(1).WriteArray(global_data[time_index])
        else:
            raise CoverageError("DefaultWriter","The given coverage is not an instance of 'TimeCoverage' or 'TimeLevelCoverage'")

    def write_variable_sea_surface_wave_stokes_drift_velocity(self):

        if (isinstance(self.coverage, TimeCoverage) or isinstance(self.coverage, TimeLevelCoverage)):

            if self.coverage.rank == 0:
                logging.info('[DefaultWriter] Gathering data from all processors')
                global_data = np.empty(
                    [self.coverage.get_t_size(type="target_global"),
                     2,
                     self.coverage.get_y_size(type="target_global"),
                     self.coverage.get_x_size(type="target_global")])
                global_data[:] = np.nan

            for time_index in range(0, self.coverage.get_t_size()):

                local_data = self.coverage.read_variable_sea_surface_wave_stokes_drift_velocity_at_time(time_index)

                if self.coverage.rank != 0:
                    self.coverage.comm.Send(np.ascontiguousarray(local_data), dest=0)
                else:
                    # Pour le proc n°1
                    global_data[self.coverage.map_mpi[self.coverage.rank]["dst_global_t"].start + time_index,
                                0,
                                self.coverage.map_mpi[self.coverage.rank]["dst_global_y"],
                                self.coverage.map_mpi[self.coverage.rank]["dst_global_x"]] = local_data[0]

                    global_data[self.coverage.map_mpi[self.coverage.rank]["dst_global_t"].start + time_index,
                                1,
                                self.coverage.map_mpi[self.coverage.rank]["dst_global_y"],
                                self.coverage.map_mpi[self.coverage.rank]["dst_global_x"]] = local_data[1]

                    # Pour les autres
                    for source in range(1, self.coverage.size):
                        recvbuf = np.empty([2, self.coverage.map_mpi[source]["dst_local_y_size"],
                                            self.coverage.map_mpi[source]["dst_local_x_size"]])
                        self.coverage.comm.Recv(recvbuf, source=source)

                        global_data[self.coverage.map_mpi[source]["dst_global_t"].start + time_index,
                                    0,
                                    self.coverage.map_mpi[source]["dst_global_y"],
                                    self.coverage.map_mpi[source]["dst_global_x"]] = recvbuf[0]

                        global_data[self.coverage.map_mpi[source]["dst_global_t"].start + time_index,
                                    1,
                                    self.coverage.map_mpi[source]["dst_global_y"],
                                    self.coverage.map_mpi[source]["dst_global_x"]] = recvbuf[1]

            self.coverage.comm.barrier()

            if self.coverage.rank == 0:

                logging.info('[DefaultWriter] Writing variable \'Surface Stokes Drift Velocity\'\'')

                file = open(self.filename, "w")

                for time_index in range(0, self.coverage.get_t_size(type="target_global")):
                    time = self.coverage.read_axis_t(type="target_global")[time_index]

                    logging.info('[DefaultWriter] Writing variable \'Surface Stokes Drift Velocity\' at time \'' + str(
                        time) + '\'')

                    for vector in range(0, 2):
                        if vector == 0:

                            file = self.driver.Create(os.path.join(self.filename, time.strftime("%Y%m%d_%H%M%S") + "_" +
                                                                   VariableDefinition.VARIABLE_NAME[
                                                                       'eastward_sea_surface_wave_stokes_drift_velocity'] + ".tiff"),
                                                      int(self.rows), int(self.cols), 1, gdal.GDT_Float64)

                            # CRS info
                            proj = osr.SpatialReference()
                            proj.SetWellKnownGeogCS("EPSG:4326")
                            file.SetProjection(proj.ExportToWkt())
                            file.SetGeoTransform(self.geotransform)
                            file.GetRasterBand(1).WriteArray(global_data[time_index][0])

                        else:

                            file = self.driver.Create(os.path.join(self.filename, time.strftime("%Y%m%d_%H%M%S") + "_" +
                                                                   VariableDefinition.VARIABLE_NAME[
                                                                       'northward_sea_surface_wave_stokes_drift_velocity'] + ".tiff"),
                                                      int(self.rows), int(self.cols), 1, gdal.GDT_Float64)

                            # CRS info
                            proj = osr.SpatialReference()
                            proj.SetWellKnownGeogCS("EPSG:4326")
                            file.SetProjection(proj.ExportToWkt())
                            file.SetGeoTransform(self.geotransform)
                            file.GetRasterBand(1).WriteArray(global_data[time_index][1])

        else:
            raise CoverageError("DefaultWriter","The given coverage is not an instance of 'TimeCoverage' or 'TimeLevelCoverage'")

    def write_variable_radiation_pressure_bernouilli_head(self):

        if (isinstance(self.coverage, TimeCoverage) or isinstance(self.coverage, TimeLevelCoverage)):

            if self.coverage.rank == 0:
                logging.info('[DefaultWriter] Gathering data from all processors')
                global_data = np.empty(
                    [self.coverage.get_t_size(type="target_global"),
                     self.coverage.get_y_size(type="target_global"),
                     self.coverage.get_x_size(type="target_global")])
                global_data[:] = np.nan

            for time_index in range(0, self.coverage.get_t_size()):

                local_data = self.coverage.read_variable_radiation_pressure_bernouilli_head_at_time(time_index)

                if self.coverage.rank != 0:
                    self.coverage.comm.Send(np.ascontiguousarray(local_data), dest=0)
                else:
                    # Pour le proc n°1
                    global_data[self.coverage.map_mpi[self.coverage.rank]["dst_global_t"].start + time_index,
                                self.coverage.map_mpi[self.coverage.rank]["dst_global_y"],
                                self.coverage.map_mpi[self.coverage.rank]["dst_global_x"]] = local_data

                    # Pour les autres
                    for source in range(1, self.coverage.size):
                        recvbuf = np.empty([self.coverage.map_mpi[source]["dst_local_y_size"],
                                            self.coverage.map_mpi[source]["dst_local_x_size"]])
                        self.coverage.comm.Recv(recvbuf, source=source)

                        global_data[self.coverage.map_mpi[source]["dst_global_t"].start + time_index,
                                    self.coverage.map_mpi[source]["dst_global_y"],
                                    self.coverage.map_mpi[source]["dst_global_x"]] = recvbuf

            self.coverage.comm.barrier()

            if self.coverage.rank == 0:

                logging.info('[DefaultWriter] Writing variable \'' + str(
                    VariableDefinition.LONG_NAME['radiation_pressure_bernouilli_head']) + '\'')

                for time_index in range(0, self.coverage.get_t_size(type="target_global")):
                    time = self.coverage.read_axis_t(type="target_global")[time_index]

                    logging.debug('[DefaultWriter] Writing variable \'' + str(VariableDefinition.LONG_NAME[
                                                                                  'radiation_pressure_bernouilli_head']) + '\' at time \'' + str(
                        time) + '\'')

                    file = self.driver.Create(os.path.join(self.filename, time.strftime("%Y%m%d_%H%M%S") + "_" +
                                                           VariableDefinition.VARIABLE_NAME[
                                                               'radiation_pressure_bernouilli_head'] + ".tiff"),
                                              int(self.rows), int(self.cols), 1, gdal.GDT_Float64)

                    # CRS info
                    proj = osr.SpatialReference()
                    proj.SetWellKnownGeogCS("EPSG:4326")
                    file.SetProjection(proj.ExportToWkt())
                    file.SetGeoTransform(self.geotransform)
                    file.GetRasterBand(1).WriteArray(global_data[time_index])
        else:
            raise CoverageError("DefaultWriter","The given coverage is not an instance of 'TimeCoverage' or 'TimeLevelCoverage'")

    def write_variable_sea_surface_wave_energy_flux_to_ocean(self):

        if (isinstance(self.coverage, TimeCoverage) or isinstance(self.coverage, TimeLevelCoverage)):

            if self.coverage.rank == 0:
                logging.info('[DefaultWriter] Gathering data from all processors')
                global_data = np.empty(
                    [self.coverage.get_t_size(type="target_global"),
                     self.coverage.get_y_size(type="target_global"),
                     self.coverage.get_x_size(type="target_global")])
                global_data[:] = np.nan

            for time_index in range(0, self.coverage.get_t_size()):

                local_data = self.coverage.read_variable_sea_surface_wave_energy_flux_to_ocean_at_time(time_index)

                if self.coverage.rank != 0:
                    self.coverage.comm.Send(np.ascontiguousarray(local_data), dest=0)
                else:
                    # Pour le proc n°1
                    global_data[self.coverage.map_mpi[self.coverage.rank]["dst_global_t"].start + time_index,
                                self.coverage.map_mpi[self.coverage.rank]["dst_global_y"],
                                self.coverage.map_mpi[self.coverage.rank]["dst_global_x"]] = local_data

                    # Pour les autres
                    for source in range(1, self.coverage.size):
                        recvbuf = np.empty([self.coverage.map_mpi[source]["dst_local_y_size"],
                                            self.coverage.map_mpi[source]["dst_local_x_size"]])
                        self.coverage.comm.Recv(recvbuf, source=source)

                        global_data[self.coverage.map_mpi[source]["dst_global_t"].start + time_index,
                                    self.coverage.map_mpi[source]["dst_global_y"],
                                    self.coverage.map_mpi[source]["dst_global_x"]] = recvbuf

            self.coverage.comm.barrier()

            if self.coverage.rank == 0:

                logging.info('[DefaultWriter] Writing variable \'' + str(
                    VariableDefinition.LONG_NAME['sea_surface_wave_energy_flux_to_ocean']) + '\'')

                for time_index in range(0, self.coverage.get_t_size(type="target_global")):
                    time = self.coverage.read_axis_t(type="target_global")[time_index]

                    logging.debug('[DefaultWriter] Writing variable \'' + str(VariableDefinition.LONG_NAME[
                                                                                  'sea_surface_wave_energy_flux_to_ocean']) + '\' at time \'' + str(
                        time) + '\'')

                    file = self.driver.Create(os.path.join(self.filename, time.strftime("%Y%m%d_%H%M%S") + "_" +
                                                           VariableDefinition.VARIABLE_NAME[
                                                               'sea_surface_wave_energy_flux_to_ocean'] + ".tiff"),
                                              int(self.rows), int(self.cols), 1, gdal.GDT_Float64)

                    # CRS info
                    proj = osr.SpatialReference()
                    proj.SetWellKnownGeogCS("EPSG:4326")
                    file.SetProjection(proj.ExportToWkt())
                    file.SetGeoTransform(self.geotransform)
                    file.GetRasterBand(1).WriteArray(global_data[time_index])
        else:
            raise CoverageError("DefaultWriter","The given coverage is not an instance of 'TimeCoverage' or 'TimeLevelCoverage'")

    #################
    # WAVES
    # Ground level
    #################
    def write_variable_sea_surface_wave_energy_dissipation_at_ground_level(self):

        if (isinstance(self.coverage, TimeCoverage) or isinstance(self.coverage, TimeLevelCoverage)):

            if self.coverage.rank == 0:
                logging.info('[DefaultWriter] Gathering data from all processors')
                global_data = np.empty(
                    [self.coverage.get_t_size(type="target_global"),
                     self.coverage.get_y_size(type="target_global"),
                     self.coverage.get_x_size(type="target_global")])
                global_data[:] = np.nan

            for time_index in range(0, self.coverage.get_t_size()):

                local_data = self.coverage.read_variable_sea_surface_wave_energy_dissipation_at_ground_level_at_time(time_index)

                if self.coverage.rank != 0:
                    self.coverage.comm.Send(np.ascontiguousarray(local_data), dest=0)
                else:
                    # Pour le proc n°1
                    global_data[self.coverage.map_mpi[self.coverage.rank]["dst_global_t"].start + time_index,
                                self.coverage.map_mpi[self.coverage.rank]["dst_global_y"],
                                self.coverage.map_mpi[self.coverage.rank]["dst_global_x"]] = local_data

                    # Pour les autres
                    for source in range(1, self.coverage.size):
                        recvbuf = np.empty([self.coverage.map_mpi[source]["dst_local_y_size"],
                                            self.coverage.map_mpi[source]["dst_local_x_size"]])
                        self.coverage.comm.Recv(recvbuf, source=source)

                        global_data[self.coverage.map_mpi[source]["dst_global_t"].start + time_index,
                                    self.coverage.map_mpi[source]["dst_global_y"],
                                    self.coverage.map_mpi[source]["dst_global_x"]] = recvbuf

            self.coverage.comm.barrier()

            if self.coverage.rank == 0:

                logging.info('[DefaultWriter] Writing variable \'' + str(
                    VariableDefinition.LONG_NAME['sea_surface_wave_energy_dissipation_at_ground_level']) + '\'')

                for time_index in range(0, self.coverage.get_t_size(type="target_global")):
                    time = self.coverage.read_axis_t(type="target_global")[time_index]

                    logging.debug('[DefaultWriter] Writing variable \'' + str(VariableDefinition.LONG_NAME[
                                                                                  'sea_surface_wave_energy_dissipation_at_ground_level']) + '\' at time \'' + str(
                        time) + '\'')

                    file = self.driver.Create(os.path.join(self.filename, time.strftime("%Y%m%d_%H%M%S") + "_" +
                                                           VariableDefinition.VARIABLE_NAME[
                                                               'sea_surface_wave_energy_dissipation_at_ground_level'] + ".tiff"),
                                              int(self.rows), int(self.cols), 1, gdal.GDT_Float64)

                    # CRS info
                    proj = osr.SpatialReference()
                    proj.SetWellKnownGeogCS("EPSG:4326")
                    file.SetProjection(proj.ExportToWkt())
                    file.SetGeoTransform(self.geotransform)
                    file.GetRasterBand(1).WriteArray(global_data[time_index])
        else:
            raise CoverageError("DefaultWriter","The given coverage is not an instance of 'TimeCoverage' or 'TimeLevelCoverage'")

    #################
    # WAVES
    # Momentum flux
    #################
    def write_variable_atmosphere_momentum_flux_to_waves(self):

        if (isinstance(self.coverage, TimeCoverage) or isinstance(self.coverage, TimeLevelCoverage)):

            if self.coverage.rank == 0:
                logging.info('[DefaultWriter] Gathering data from all processors')
                global_data = np.empty(
                    [self.coverage.get_t_size(type="target_global"),
                     2,
                     self.coverage.get_y_size(type="target_global"),
                     self.coverage.get_x_size(type="target_global")])
                global_data[:] = np.nan

            for time_index in range(0, self.coverage.get_t_size()):

                local_data = self.coverage.read_variable_atmosphere_momentum_flux_to_waves_at_time(time_index)

                if self.coverage.rank != 0:
                    self.coverage.comm.Send(np.ascontiguousarray(local_data), dest=0)
                else:
                    # Pour le proc n°1
                    global_data[self.coverage.map_mpi[self.coverage.rank]["dst_global_t"].start + time_index,
                                0,
                                self.coverage.map_mpi[self.coverage.rank]["dst_global_y"],
                                self.coverage.map_mpi[self.coverage.rank]["dst_global_x"]] = local_data[0]

                    global_data[self.coverage.map_mpi[self.coverage.rank]["dst_global_t"].start + time_index,
                                1,
                                self.coverage.map_mpi[self.coverage.rank]["dst_global_y"],
                                self.coverage.map_mpi[self.coverage.rank]["dst_global_x"]] = local_data[1]

                    # Pour les autres
                    for source in range(1, self.coverage.size):
                        recvbuf = np.empty([2, self.coverage.map_mpi[source]["dst_local_y_size"],
                                            self.coverage.map_mpi[source]["dst_local_x_size"]])
                        self.coverage.comm.Recv(recvbuf, source=source)

                        global_data[self.coverage.map_mpi[source]["dst_global_t"].start + time_index,
                                    0,
                                    self.coverage.map_mpi[source]["dst_global_y"],
                                    self.coverage.map_mpi[source]["dst_global_x"]] = recvbuf[0]

                        global_data[self.coverage.map_mpi[source]["dst_global_t"].start + time_index,
                                    1,
                                    self.coverage.map_mpi[source]["dst_global_y"],
                                    self.coverage.map_mpi[source]["dst_global_x"]] = recvbuf[1]

            self.coverage.comm.barrier()

            if self.coverage.rank == 0:

                logging.info('[DefaultWriter] Writing variable \'Atmosphere Momentum Flux to Waves\'\'')

                file = open(self.filename, "w")

                for time_index in range(0, self.coverage.get_t_size(type="target_global")):
                    time = self.coverage.read_axis_t(type="target_global")[time_index]

                    logging.info(
                        '[DefaultWriter] Writing variable \'Atmosphere Momentum Flux to Waves\' at time \'' + str(
                            time) + '\'')

                    for vector in range(0, 2):
                        if vector == 0:

                            file = self.driver.Create(os.path.join(self.filename, time.strftime("%Y%m%d_%H%M%S") + "_" +
                                                                   VariableDefinition.VARIABLE_NAME[
                                                                       'eastward_atmosphere_momentum_flux_to_waves'] + ".tiff"),
                                                      int(self.rows), int(self.cols), 1, gdal.GDT_Float64)

                            # CRS info
                            proj = osr.SpatialReference()
                            proj.SetWellKnownGeogCS("EPSG:4326")
                            file.SetProjection(proj.ExportToWkt())
                            file.SetGeoTransform(self.geotransform)
                            file.GetRasterBand(1).WriteArray(global_data[time_index][0])

                        else:

                            file = self.driver.Create(os.path.join(self.filename, time.strftime("%Y%m%d_%H%M%S") + "_" +
                                                                   VariableDefinition.VARIABLE_NAME[
                                                                       'northward_atmosphere_momentum_flux_to_waves'] + ".tiff"),
                                                      int(self.rows), int(self.cols), 1, gdal.GDT_Float64)

                            # CRS info
                            proj = osr.SpatialReference()
                            proj.SetWellKnownGeogCS("EPSG:4326")
                            file.SetProjection(proj.ExportToWkt())
                            file.SetGeoTransform(self.geotransform)
                            file.GetRasterBand(1).WriteArray(global_data[time_index][1])

        else:
            raise CoverageError("DefaultWriter","The given coverage is not an instance of 'TimeCoverage' or 'TimeLevelCoverage'")

    def write_variable_waves_momentum_flux_to_ocean(self):

        if (isinstance(self.coverage, TimeCoverage) or isinstance(self.coverage, TimeLevelCoverage)):

            if self.coverage.rank == 0:
                logging.info('[DefaultWriter] Gathering data from all processors')
                global_data = np.empty(
                    [self.coverage.get_t_size(type="target_global"),
                     2,
                     self.coverage.get_y_size(type="target_global"),
                     self.coverage.get_x_size(type="target_global")])
                global_data[:] = np.nan

            for time_index in range(0, self.coverage.get_t_size()):

                local_data = self.coverage.read_variable_waves_momentum_flux_to_ocean_at_time(time_index)

                if self.coverage.rank != 0:
                    self.coverage.comm.Send(np.ascontiguousarray(local_data), dest=0)
                else:
                    # Pour le proc n°1
                    global_data[self.coverage.map_mpi[self.coverage.rank]["dst_global_t"].start + time_index,
                                0,
                                self.coverage.map_mpi[self.coverage.rank]["dst_global_y"],
                                self.coverage.map_mpi[self.coverage.rank]["dst_global_x"]] = local_data[0]

                    global_data[self.coverage.map_mpi[self.coverage.rank]["dst_global_t"].start + time_index,
                                1,
                                self.coverage.map_mpi[self.coverage.rank]["dst_global_y"],
                                self.coverage.map_mpi[self.coverage.rank]["dst_global_x"]] = local_data[1]

                    # Pour les autres
                    for source in range(1, self.coverage.size):
                        recvbuf = np.empty([2, self.coverage.map_mpi[source]["dst_local_y_size"],
                                            self.coverage.map_mpi[source]["dst_local_x_size"]])
                        self.coverage.comm.Recv(recvbuf, source=source)

                        global_data[self.coverage.map_mpi[source]["dst_global_t"].start + time_index,
                                    0,
                                    self.coverage.map_mpi[source]["dst_global_y"],
                                    self.coverage.map_mpi[source]["dst_global_x"]] = recvbuf[0]

                        global_data[self.coverage.map_mpi[source]["dst_global_t"].start + time_index,
                                    1,
                                    self.coverage.map_mpi[source]["dst_global_y"],
                                    self.coverage.map_mpi[source]["dst_global_x"]] = recvbuf[1]

            self.coverage.comm.barrier()

            if self.coverage.rank == 0:

                logging.info('[DefaultWriter] Writing variable \'Waves Momentum Flux To Ocean\'\'')

                file = open(self.filename, "w")

                for time_index in range(0, self.coverage.get_t_size(type="target_global")):
                    time = self.coverage.read_axis_t(type="target_global")[time_index]

                    logging.info('[DefaultWriter] Writing variable \'Waves Momentum Flux To Ocean\' at time \''+str(time)+'\'')

                    for vector in range(0, 2):
                        if vector == 0:

                            file = self.driver.Create(os.path.join(self.filename, time.strftime("%Y%m%d_%H%M%S") + "_" +
                                                                   VariableDefinition.VARIABLE_NAME[
                                                                       'eastward_waves_momentum_flux_to_ocean'] + ".tiff"),
                                                      int(self.rows), int(self.cols), 1, gdal.GDT_Float64)

                            # CRS info
                            proj = osr.SpatialReference()
                            proj.SetWellKnownGeogCS("EPSG:4326")
                            file.SetProjection(proj.ExportToWkt())
                            file.SetGeoTransform(self.geotransform)
                            file.GetRasterBand(1).WriteArray(global_data[time_index][0])

                        else:

                            file = self.driver.Create(os.path.join(self.filename, time.strftime("%Y%m%d_%H%M%S") + "_" +
                                                                   VariableDefinition.VARIABLE_NAME[
                                                                       'northward_waves_momentum_flux_to_ocean'] + ".tiff"),
                                                      int(self.rows), int(self.cols), 1, gdal.GDT_Float64)

                            # CRS info
                            proj = osr.SpatialReference()
                            proj.SetWellKnownGeogCS("EPSG:4326")
                            file.SetProjection(proj.ExportToWkt())
                            file.SetGeoTransform(self.geotransform)
                            file.GetRasterBand(1).WriteArray(global_data[time_index][1])

        else:
            raise CoverageError("DefaultWriter","The given coverage is not an instance of 'TimeCoverage' or 'TimeLevelCoverage'")

    #################
    # METEO
    # 2D
    #################

    #################
    # METEO
    # Surface air
    #################
    def write_variable_sea_surface_air_pressure(self):

        if (isinstance(self.coverage, TimeCoverage) or isinstance(self.coverage, TimeLevelCoverage)):

            if self.coverage.rank == 0:
                logging.info('[DefaultWriter] Gathering data from all processors')
                global_data = np.empty(
                    [self.coverage.get_t_size(type="target_global"),
                     self.coverage.get_y_size(type="target_global"),
                     self.coverage.get_x_size(type="target_global")])
                global_data[:] = np.nan

            for time_index in range(0, self.coverage.get_t_size()):

                local_data = self.coverage.read_variable_sea_surface_air_pressure_at_time(
                    time_index)

                if self.coverage.rank != 0:
                    self.coverage.comm.Send(np.ascontiguousarray(local_data), dest=0)
                else:
                    # Pour le proc n°1
                    global_data[self.coverage.map_mpi[self.coverage.rank]["dst_global_t"].start + time_index,
                                self.coverage.map_mpi[self.coverage.rank]["dst_global_y"],
                                self.coverage.map_mpi[self.coverage.rank]["dst_global_x"]] = local_data

                    # Pour les autres
                    for source in range(1, self.coverage.size):
                        recvbuf = np.empty([self.coverage.map_mpi[source]["dst_local_y_size"],
                                            self.coverage.map_mpi[source]["dst_local_x_size"]])
                        self.coverage.comm.Recv(recvbuf, source=source)

                        global_data[self.coverage.map_mpi[source]["dst_global_t"].start + time_index,
                                    self.coverage.map_mpi[source]["dst_global_y"],
                                    self.coverage.map_mpi[source]["dst_global_x"]] = recvbuf

            self.coverage.comm.barrier()

            if self.coverage.rank == 0:

                logging.info('[DefaultWriter] Writing variable \'' + str(
                    VariableDefinition.LONG_NAME['sea_surface_air_pressure']) + '\'')

                for time_index in range(0, self.coverage.get_t_size(type="target_global")):
                    time = self.coverage.read_axis_t(type="target_global")[time_index]

                    logging.debug('[DefaultWriter] Writing variable \'' + str(VariableDefinition.LONG_NAME[
                                                                                  'sea_surface_air_pressure']) + '\' at time \'' + str(
                        time) + '\'')

                    file = self.driver.Create(os.path.join(self.filename, time.strftime("%Y%m%d_%H%M%S") + "_" +
                                                           VariableDefinition.VARIABLE_NAME[
                                                               'sea_surface_air_pressure'] + ".tiff"),
                                              int(self.rows), int(self.cols), 1, gdal.GDT_Float64)

                    # CRS info
                    proj = osr.SpatialReference()
                    proj.SetWellKnownGeogCS("EPSG:4326")
                    file.SetProjection(proj.ExportToWkt())
                    file.SetGeoTransform(self.geotransform)
                    file.GetRasterBand(1).WriteArray(global_data[time_index])
        else:
            raise CoverageError("DefaultWriter",
                                "The given coverage is not an instance of 'TimeCoverage' or 'TimeLevelCoverage'")

    def write_variable_surface_air_temperature(self):

        if (isinstance(self.coverage, TimeCoverage) or isinstance(self.coverage, TimeLevelCoverage)):

            if self.coverage.rank == 0:
                logging.info('[DefaultWriter] Gathering data from all processors')
                global_data = np.empty(
                    [self.coverage.get_t_size(type="target_global"),
                     self.coverage.get_y_size(type="target_global"),
                     self.coverage.get_x_size(type="target_global")])
                global_data[:] = np.nan

            for time_index in range(0, self.coverage.get_t_size()):

                local_data = self.coverage.read_variable_surface_air_temperature_at_time(
                    time_index)

                if self.coverage.rank != 0:
                    self.coverage.comm.Send(np.ascontiguousarray(local_data), dest=0)
                else:
                    # Pour le proc n°1
                    global_data[self.coverage.map_mpi[self.coverage.rank]["dst_global_t"].start + time_index,
                                self.coverage.map_mpi[self.coverage.rank]["dst_global_y"],
                                self.coverage.map_mpi[self.coverage.rank]["dst_global_x"]] = local_data

                    # Pour les autres
                    for source in range(1, self.coverage.size):
                        recvbuf = np.empty([self.coverage.map_mpi[source]["dst_local_y_size"],
                                            self.coverage.map_mpi[source]["dst_local_x_size"]])
                        self.coverage.comm.Recv(recvbuf, source=source)

                        global_data[self.coverage.map_mpi[source]["dst_global_t"].start + time_index,
                                    self.coverage.map_mpi[source]["dst_global_y"],
                                    self.coverage.map_mpi[source]["dst_global_x"]] = recvbuf

            self.coverage.comm.barrier()

            if self.coverage.rank == 0:

                logging.info('[DefaultWriter] Writing variable \'' + str(
                    VariableDefinition.LONG_NAME['surface_air_temperature']) + '\'')

                for time_index in range(0, self.coverage.get_t_size(type="target_global")):
                    time = self.coverage.read_axis_t(type="target_global")[time_index]

                    logging.debug('[DefaultWriter] Writing variable \'' + str(VariableDefinition.LONG_NAME[
                                                                                  'surface_air_temperature']) + '\' at time \'' + str(
                        time) + '\'')

                    file = self.driver.Create(os.path.join(self.filename, time.strftime("%Y%m%d_%H%M%S") + "_" +
                                                           VariableDefinition.VARIABLE_NAME[
                                                               'surface_air_temperature'] + ".tiff"),
                                              int(self.rows), int(self.cols), 1, gdal.GDT_Float64)

                    # CRS info
                    proj = osr.SpatialReference()
                    proj.SetWellKnownGeogCS("EPSG:4326")
                    file.SetProjection(proj.ExportToWkt())
                    file.SetGeoTransform(self.geotransform)
                    file.GetRasterBand(1).WriteArray(global_data[time_index])
        else:
            raise CoverageError("DefaultWriter",
                                "The given coverage is not an instance of 'TimeCoverage' or 'TimeLevelCoverage'")

    def write_variable_dew_point_temperature(self):

        if (isinstance(self.coverage, TimeCoverage) or isinstance(self.coverage, TimeLevelCoverage)):

            if self.coverage.rank == 0:
                logging.info('[DefaultWriter] Gathering data from all processors')
                global_data = np.empty(
                    [self.coverage.get_t_size(type="target_global"),
                     self.coverage.get_y_size(type="target_global"),
                     self.coverage.get_x_size(type="target_global")])
                global_data[:] = np.nan

            for time_index in range(0, self.coverage.get_t_size()):

                local_data = self.coverage.read_variable_dew_point_temperature_at_time(
                    time_index)

                if self.coverage.rank != 0:
                    self.coverage.comm.Send(np.ascontiguousarray(local_data), dest=0)
                else:
                    # Pour le proc n°1
                    global_data[self.coverage.map_mpi[self.coverage.rank]["dst_global_t"].start + time_index,
                                self.coverage.map_mpi[self.coverage.rank]["dst_global_y"],
                                self.coverage.map_mpi[self.coverage.rank]["dst_global_x"]] = local_data

                    # Pour les autres
                    for source in range(1, self.coverage.size):
                        recvbuf = np.empty([self.coverage.map_mpi[source]["dst_local_y_size"],
                                            self.coverage.map_mpi[source]["dst_local_x_size"]])
                        self.coverage.comm.Recv(recvbuf, source=source)

                        global_data[self.coverage.map_mpi[source]["dst_global_t"].start + time_index,
                                    self.coverage.map_mpi[source]["dst_global_y"],
                                    self.coverage.map_mpi[source]["dst_global_x"]] = recvbuf

            self.coverage.comm.barrier()

            if self.coverage.rank == 0:

                logging.info('[DefaultWriter] Writing variable \'' + str(
                    VariableDefinition.LONG_NAME['dew_point_temperature']) + '\'')

                for time_index in range(0, self.coverage.get_t_size(type="target_global")):
                    time = self.coverage.read_axis_t(type="target_global")[time_index]

                    logging.debug('[DefaultWriter] Writing variable \'' + str(VariableDefinition.LONG_NAME[
                                                                                  'dew_point_temperature']) + '\' at time \'' + str(
                        time) + '\'')

                    file = self.driver.Create(os.path.join(self.filename, time.strftime("%Y%m%d_%H%M%S") + "_" +
                                                           VariableDefinition.VARIABLE_NAME[
                                                               'dew_point_temperature'] + ".tiff"),
                                              int(self.rows), int(self.cols), 1, gdal.GDT_Float64)

                    # CRS info
                    proj = osr.SpatialReference()
                    proj.SetWellKnownGeogCS("EPSG:4326")
                    file.SetProjection(proj.ExportToWkt())
                    file.SetGeoTransform(self.geotransform)
                    file.GetRasterBand(1).WriteArray(global_data[time_index])
        else:
            raise CoverageError("DefaultWriter",
                                "The given coverage is not an instance of 'TimeCoverage' or 'TimeLevelCoverage'")

    def write_variable_rainfall_amount(self):

        if (isinstance(self.coverage, TimeCoverage) or isinstance(self.coverage, TimeLevelCoverage)):

            if self.coverage.rank == 0:
                logging.info('[DefaultWriter] Gathering data from all processors')
                global_data = np.empty(
                    [self.coverage.get_t_size(type="target_global"),
                     self.coverage.get_y_size(type="target_global"),
                     self.coverage.get_x_size(type="target_global")])
                global_data[:] = np.nan

            for time_index in range(0, self.coverage.get_t_size()):

                local_data = self.coverage.read_variable_rainfall_amount_at_time(
                    time_index)

                if self.coverage.rank != 0:
                    self.coverage.comm.Send(np.ascontiguousarray(local_data), dest=0)
                else:
                    # Pour le proc n°1
                    global_data[self.coverage.map_mpi[self.coverage.rank]["dst_global_t"].start + time_index,
                                self.coverage.map_mpi[self.coverage.rank]["dst_global_y"],
                                self.coverage.map_mpi[self.coverage.rank]["dst_global_x"]] = local_data

                    # Pour les autres
                    for source in range(1, self.coverage.size):
                        recvbuf = np.empty([self.coverage.map_mpi[source]["dst_local_y_size"],
                                            self.coverage.map_mpi[source]["dst_local_x_size"]])
                        self.coverage.comm.Recv(recvbuf, source=source)

                        global_data[self.coverage.map_mpi[source]["dst_global_t"].start + time_index,
                                    self.coverage.map_mpi[source]["dst_global_y"],
                                    self.coverage.map_mpi[source]["dst_global_x"]] = recvbuf

            self.coverage.comm.barrier()

            if self.coverage.rank == 0:

                logging.info('[DefaultWriter] Writing variable \'' + str(
                    VariableDefinition.LONG_NAME['rainfall_amount']) + '\'')

                for time_index in range(0, self.coverage.get_t_size(type="target_global")):
                    time = self.coverage.read_axis_t(type="target_global")[time_index]

                    logging.debug('[DefaultWriter] Writing variable \'' + str(VariableDefinition.LONG_NAME[
                                                                                  'rainfall_amount']) + '\' at time \'' + str(
                        time) + '\'')

                    file = self.driver.Create(os.path.join(self.filename, time.strftime("%Y%m%d_%H%M%S") + "_" +
                                                           VariableDefinition.VARIABLE_NAME[
                                                               'rainfall_amount'] + ".tiff"),
                                              int(self.rows), int(self.cols), 1, gdal.GDT_Float64)

                    # CRS info
                    proj = osr.SpatialReference()
                    proj.SetWellKnownGeogCS("EPSG:4326")
                    file.SetProjection(proj.ExportToWkt())
                    file.SetGeoTransform(self.geotransform)
                    file.GetRasterBand(1).WriteArray(global_data[time_index])
        else:
            raise CoverageError("DefaultWriter",
                                "The given coverage is not an instance of 'TimeCoverage' or 'TimeLevelCoverage'")

    def write_variable_surface_downward_sensible_heat_flux(self):

        if (isinstance(self.coverage, TimeCoverage) or isinstance(self.coverage, TimeLevelCoverage)):

            if self.coverage.rank == 0:
                logging.info('[DefaultWriter] Gathering data from all processors')
                global_data = np.empty(
                    [self.coverage.get_t_size(type="target_global"),
                     self.coverage.get_y_size(type="target_global"),
                     self.coverage.get_x_size(type="target_global")])
                global_data[:] = np.nan

            for time_index in range(0, self.coverage.get_t_size()):

                local_data = self.coverage.read_variable_surface_downward_sensible_heat_flux_at_time(
                    time_index)

                if self.coverage.rank != 0:
                    self.coverage.comm.Send(np.ascontiguousarray(local_data), dest=0)
                else:
                    # Pour le proc n°1
                    global_data[self.coverage.map_mpi[self.coverage.rank]["dst_global_t"].start + time_index,
                                self.coverage.map_mpi[self.coverage.rank]["dst_global_y"],
                                self.coverage.map_mpi[self.coverage.rank]["dst_global_x"]] = local_data

                    # Pour les autres
                    for source in range(1, self.coverage.size):
                        recvbuf = np.empty([self.coverage.map_mpi[source]["dst_local_y_size"],
                                            self.coverage.map_mpi[source]["dst_local_x_size"]])
                        self.coverage.comm.Recv(recvbuf, source=source)

                        global_data[self.coverage.map_mpi[source]["dst_global_t"].start + time_index,
                                    self.coverage.map_mpi[source]["dst_global_y"],
                                    self.coverage.map_mpi[source]["dst_global_x"]] = recvbuf

            self.coverage.comm.barrier()

            if self.coverage.rank == 0:

                logging.info('[DefaultWriter] Writing variable \'' + str(
                    VariableDefinition.LONG_NAME['surface_downward_sensible_heat_flux']) + '\'')

                for time_index in range(0, self.coverage.get_t_size(type="target_global")):
                    time = self.coverage.read_axis_t(type="target_global")[time_index]

                    logging.debug('[DefaultWriter] Writing variable \'' + str(VariableDefinition.LONG_NAME[
                                                                                  'surface_downward_sensible_heat_flux']) + '\' at time \'' + str(
                        time) + '\'')

                    file = self.driver.Create(os.path.join(self.filename, time.strftime("%Y%m%d_%H%M%S") + "_" +
                                                           VariableDefinition.VARIABLE_NAME[
                                                               'surface_downward_sensible_heat_flux'] + ".tiff"),
                                              int(self.rows), int(self.cols), 1, gdal.GDT_Float64)

                    # CRS info
                    proj = osr.SpatialReference()
                    proj.SetWellKnownGeogCS("EPSG:4326")
                    file.SetProjection(proj.ExportToWkt())
                    file.SetGeoTransform(self.geotransform)
                    file.GetRasterBand(1).WriteArray(global_data[time_index])
        else:
            raise CoverageError("DefaultWriter",
                                "The given coverage is not an instance of 'TimeCoverage' or 'TimeLevelCoverage'")

    def write_variable_surface_downward_latent_heat_flux(self):

        if (isinstance(self.coverage, TimeCoverage) or isinstance(self.coverage, TimeLevelCoverage)):

            if self.coverage.rank == 0:
                logging.info('[DefaultWriter] Gathering data from all processors')
                global_data = np.empty(
                    [self.coverage.get_t_size(type="target_global"),
                     self.coverage.get_y_size(type="target_global"),
                     self.coverage.get_x_size(type="target_global")])
                global_data[:] = np.nan

            for time_index in range(0, self.coverage.get_t_size()):

                local_data = self.coverage.read_variable_surface_downward_latent_heat_flux_at_time(
                    time_index)

                if self.coverage.rank != 0:
                    self.coverage.comm.Send(np.ascontiguousarray(local_data), dest=0)
                else:
                    # Pour le proc n°1
                    global_data[self.coverage.map_mpi[self.coverage.rank]["dst_global_t"].start + time_index,
                                self.coverage.map_mpi[self.coverage.rank]["dst_global_y"],
                                self.coverage.map_mpi[self.coverage.rank]["dst_global_x"]] = local_data

                    # Pour les autres
                    for source in range(1, self.coverage.size):
                        recvbuf = np.empty([self.coverage.map_mpi[source]["dst_local_y_size"],
                                            self.coverage.map_mpi[source]["dst_local_x_size"]])
                        self.coverage.comm.Recv(recvbuf, source=source)

                        global_data[self.coverage.map_mpi[source]["dst_global_t"].start + time_index,
                                    self.coverage.map_mpi[source]["dst_global_y"],
                                    self.coverage.map_mpi[source]["dst_global_x"]] = recvbuf

            self.coverage.comm.barrier()

            if self.coverage.rank == 0:

                logging.info('[DefaultWriter] Writing variable \'' + str(
                    VariableDefinition.LONG_NAME['surface_downward_latent_heat_flux']) + '\'')

                for time_index in range(0, self.coverage.get_t_size(type="target_global")):
                    time = self.coverage.read_axis_t(type="target_global")[time_index]

                    logging.debug('[DefaultWriter] Writing variable \'' + str(VariableDefinition.LONG_NAME[
                                                                                  'surface_downward_latent_heat_flux']) + '\' at time \'' + str(
                        time) + '\'')

                    file = self.driver.Create(os.path.join(self.filename, time.strftime("%Y%m%d_%H%M%S") + "_" +
                                                           VariableDefinition.VARIABLE_NAME[
                                                               'surface_downward_latent_heat_flux'] + ".tiff"),
                                              int(self.rows), int(self.cols), 1, gdal.GDT_Float64)

                    # CRS info
                    proj = osr.SpatialReference()
                    proj.SetWellKnownGeogCS("EPSG:4326")
                    file.SetProjection(proj.ExportToWkt())
                    file.SetGeoTransform(self.geotransform)
                    file.GetRasterBand(1).WriteArray(global_data[time_index])
        else:
            raise CoverageError("DefaultWriter",
                                "The given coverage is not an instance of 'TimeCoverage' or 'TimeLevelCoverage'")

    def write_variable_surface_downward_solar_radiation(self):

        if (isinstance(self.coverage, TimeCoverage) or isinstance(self.coverage, TimeLevelCoverage)):

            if self.coverage.rank == 0:
                logging.info('[DefaultWriter] Gathering data from all processors')
                global_data = np.empty(
                    [self.coverage.get_t_size(type="target_global"),
                     self.coverage.get_y_size(type="target_global"),
                     self.coverage.get_x_size(type="target_global")])
                global_data[:] = np.nan

            for time_index in range(0, self.coverage.get_t_size()):

                local_data = self.coverage.read_variable_surface_downward_solar_radiation_at_time(
                    time_index)

                if self.coverage.rank != 0:
                    self.coverage.comm.Send(np.ascontiguousarray(local_data), dest=0)
                else:
                    # Pour le proc n°1
                    global_data[self.coverage.map_mpi[self.coverage.rank]["dst_global_t"].start + time_index,
                                self.coverage.map_mpi[self.coverage.rank]["dst_global_y"],
                                self.coverage.map_mpi[self.coverage.rank]["dst_global_x"]] = local_data

                    # Pour les autres
                    for source in range(1, self.coverage.size):
                        recvbuf = np.empty([self.coverage.map_mpi[source]["dst_local_y_size"],
                                            self.coverage.map_mpi[source]["dst_local_x_size"]])
                        self.coverage.comm.Recv(recvbuf, source=source)

                        global_data[self.coverage.map_mpi[source]["dst_global_t"].start + time_index,
                                    self.coverage.map_mpi[source]["dst_global_y"],
                                    self.coverage.map_mpi[source]["dst_global_x"]] = recvbuf

            self.coverage.comm.barrier()

            if self.coverage.rank == 0:

                logging.info('[DefaultWriter] Writing variable \'' + str(
                    VariableDefinition.LONG_NAME['surface_downwards_solar_radiation']) + '\'')

                for time_index in range(0, self.coverage.get_t_size(type="target_global")):
                    time = self.coverage.read_axis_t(type="target_global")[time_index]

                    logging.debug('[DefaultWriter] Writing variable \'' + str(VariableDefinition.LONG_NAME[
                                                                                  'surface_downwards_solar_radiation']) + '\' at time \'' + str(
                        time) + '\'')

                    file = self.driver.Create(os.path.join(self.filename, time.strftime("%Y%m%d_%H%M%S") + "_" +
                                                           VariableDefinition.VARIABLE_NAME[
                                                               'surface_downwards_solar_radiation'] + ".tiff"),
                                              int(self.rows), int(self.cols), 1, gdal.GDT_Float64)

                    # CRS info
                    proj = osr.SpatialReference()
                    proj.SetWellKnownGeogCS("EPSG:4326")
                    file.SetProjection(proj.ExportToWkt())
                    file.SetGeoTransform(self.geotransform)
                    file.GetRasterBand(1).WriteArray(global_data[time_index])
        else:
            raise CoverageError("DefaultWriter",
                                "The given coverage is not an instance of 'TimeCoverage' or 'TimeLevelCoverage'")

    def write_variable_surface_downward_thermal_radiation(self):

        if (isinstance(self.coverage, TimeCoverage) or isinstance(self.coverage, TimeLevelCoverage)):

            if self.coverage.rank == 0:
                logging.info('[DefaultWriter] Gathering data from all processors')
                global_data = np.empty(
                    [self.coverage.get_t_size(type="target_global"),
                     self.coverage.get_y_size(type="target_global"),
                     self.coverage.get_x_size(type="target_global")])
                global_data[:] = np.nan

            for time_index in range(0, self.coverage.get_t_size()):

                local_data = self.coverage.read_variable_surface_downward_thermal_radiation_at_time(
                    time_index)

                if self.coverage.rank != 0:
                    self.coverage.comm.Send(np.ascontiguousarray(local_data), dest=0)
                else:
                    # Pour le proc n°1
                    global_data[self.coverage.map_mpi[self.coverage.rank]["dst_global_t"].start + time_index,
                                self.coverage.map_mpi[self.coverage.rank]["dst_global_y"],
                                self.coverage.map_mpi[self.coverage.rank]["dst_global_x"]] = local_data

                    # Pour les autres
                    for source in range(1, self.coverage.size):
                        recvbuf = np.empty([self.coverage.map_mpi[source]["dst_local_y_size"],
                                            self.coverage.map_mpi[source]["dst_local_x_size"]])
                        self.coverage.comm.Recv(recvbuf, source=source)

                        global_data[self.coverage.map_mpi[source]["dst_global_t"].start + time_index,
                                    self.coverage.map_mpi[source]["dst_global_y"],
                                    self.coverage.map_mpi[source]["dst_global_x"]] = recvbuf

            self.coverage.comm.barrier()

            if self.coverage.rank == 0:

                logging.info('[DefaultWriter] Writing variable \'' + str(
                    VariableDefinition.LONG_NAME['surface_downwards_thermal_radiation']) + '\'')

                for time_index in range(0, self.coverage.get_t_size(type="target_global")):
                    time = self.coverage.read_axis_t(type="target_global")[time_index]

                    logging.debug('[DefaultWriter] Writing variable \'' + str(VariableDefinition.LONG_NAME[
                                                                                  'surface_downwards_thermal_radiation']) + '\' at time \'' + str(
                        time) + '\'')

                    file = self.driver.Create(os.path.join(self.filename, time.strftime("%Y%m%d_%H%M%S") + "_" +
                                                           VariableDefinition.VARIABLE_NAME[
                                                               'surface_downwards_thermal_radiation'] + ".tiff"),
                                              int(self.rows), int(self.cols), 1, gdal.GDT_Float64)

                    # CRS info
                    proj = osr.SpatialReference()
                    proj.SetWellKnownGeogCS("EPSG:4326")
                    file.SetProjection(proj.ExportToWkt())
                    file.SetGeoTransform(self.geotransform)
                    file.GetRasterBand(1).WriteArray(global_data[time_index])
        else:
            raise CoverageError("DefaultWriter",
                                "The given coverage is not an instance of 'TimeCoverage' or 'TimeLevelCoverage'")

    def write_variable_surface_solar_radiation(self):

        if (isinstance(self.coverage, TimeCoverage) or isinstance(self.coverage, TimeLevelCoverage)):

            if self.coverage.rank == 0:
                logging.info('[DefaultWriter] Gathering data from all processors')
                global_data = np.empty(
                    [self.coverage.get_t_size(type="target_global"),
                     self.coverage.get_y_size(type="target_global"),
                     self.coverage.get_x_size(type="target_global")])
                global_data[:] = np.nan

            for time_index in range(0, self.coverage.get_t_size()):

                local_data = self.coverage.read_variable_surface_solar_radiation_at_time(
                    time_index)

                if self.coverage.rank != 0:
                    self.coverage.comm.Send(np.ascontiguousarray(local_data), dest=0)
                else:
                    # Pour le proc n°1
                    global_data[self.coverage.map_mpi[self.coverage.rank]["dst_global_t"].start + time_index,
                                self.coverage.map_mpi[self.coverage.rank]["dst_global_y"],
                                self.coverage.map_mpi[self.coverage.rank]["dst_global_x"]] = local_data

                    # Pour les autres
                    for source in range(1, self.coverage.size):
                        recvbuf = np.empty([self.coverage.map_mpi[source]["dst_local_y_size"],
                                            self.coverage.map_mpi[source]["dst_local_x_size"]])
                        self.coverage.comm.Recv(recvbuf, source=source)

                        global_data[self.coverage.map_mpi[source]["dst_global_t"].start + time_index,
                                    self.coverage.map_mpi[source]["dst_global_y"],
                                    self.coverage.map_mpi[source]["dst_global_x"]] = recvbuf

            self.coverage.comm.barrier()

            if self.coverage.rank == 0:

                logging.info('[DefaultWriter] Writing variable \'' + str(
                    VariableDefinition.LONG_NAME['surface_solar_radiation']) + '\'')

                for time_index in range(0, self.coverage.get_t_size(type="target_global")):
                    time = self.coverage.read_axis_t(type="target_global")[time_index]

                    logging.debug('[DefaultWriter] Writing variable \'' + str(VariableDefinition.LONG_NAME[
                                                                                  'surface_solar_radiation']) + '\' at time \'' + str(
                        time) + '\'')

                    file = self.driver.Create(os.path.join(self.filename, time.strftime("%Y%m%d_%H%M%S") + "_" +
                                                           VariableDefinition.VARIABLE_NAME[
                                                               'surface_solar_radiation'] + ".tiff"),
                                              int(self.rows), int(self.cols), 1, gdal.GDT_Float64)

                    # CRS info
                    proj = osr.SpatialReference()
                    proj.SetWellKnownGeogCS("EPSG:4326")
                    file.SetProjection(proj.ExportToWkt())
                    file.SetGeoTransform(self.geotransform)
                    file.GetRasterBand(1).WriteArray(global_data[time_index])
        else:
            raise CoverageError("DefaultWriter",
                                "The given coverage is not an instance of 'TimeCoverage' or 'TimeLevelCoverage'")

    def write_variable_surface_thermal_radiation(self):

        if (isinstance(self.coverage, TimeCoverage) or isinstance(self.coverage, TimeLevelCoverage)):

            if self.coverage.rank == 0:
                logging.info('[DefaultWriter] Gathering data from all processors')
                global_data = np.empty(
                    [self.coverage.get_t_size(type="target_global"),
                     self.coverage.get_y_size(type="target_global"),
                     self.coverage.get_x_size(type="target_global")])
                global_data[:] = np.nan

            for time_index in range(0, self.coverage.get_t_size()):

                local_data = self.coverage.read_variable_surface_thermal_radiation_at_time(
                    time_index)

                if self.coverage.rank != 0:
                    self.coverage.comm.Send(np.ascontiguousarray(local_data), dest=0)
                else:
                    # Pour le proc n°1
                    global_data[self.coverage.map_mpi[self.coverage.rank]["dst_global_t"].start + time_index,
                                self.coverage.map_mpi[self.coverage.rank]["dst_global_y"],
                                self.coverage.map_mpi[self.coverage.rank]["dst_global_x"]] = local_data

                    # Pour les autres
                    for source in range(1, self.coverage.size):
                        recvbuf = np.empty([self.coverage.map_mpi[source]["dst_local_y_size"],
                                            self.coverage.map_mpi[source]["dst_local_x_size"]])
                        self.coverage.comm.Recv(recvbuf, source=source)

                        global_data[self.coverage.map_mpi[source]["dst_global_t"].start + time_index,
                                    self.coverage.map_mpi[source]["dst_global_y"],
                                    self.coverage.map_mpi[source]["dst_global_x"]] = recvbuf

            self.coverage.comm.barrier()

            if self.coverage.rank == 0:

                logging.info('[DefaultWriter] Writing variable \'' + str(
                    VariableDefinition.LONG_NAME['surface_thermal_radiation']) + '\'')

                for time_index in range(0, self.coverage.get_t_size(type="target_global")):
                    time = self.coverage.read_axis_t(type="target_global")[time_index]

                    logging.debug('[DefaultWriter] Writing variable \'' + str(VariableDefinition.LONG_NAME[
                                                                                  'surface_thermal_radiation']) + '\' at time \'' + str(
                        time) + '\'')

                    file = self.driver.Create(os.path.join(self.filename, time.strftime("%Y%m%d_%H%M%S") + "_" +
                                                           VariableDefinition.VARIABLE_NAME[
                                                               'surface_thermal_radiation'] + ".tiff"),
                                              int(self.rows), int(self.cols), 1, gdal.GDT_Float64)

                    # CRS info
                    proj = osr.SpatialReference()
                    proj.SetWellKnownGeogCS("EPSG:4326")
                    file.SetProjection(proj.ExportToWkt())
                    file.SetGeoTransform(self.geotransform)
                    file.GetRasterBand(1).WriteArray(global_data[time_index])
        else:
            raise CoverageError("DefaultWriter",
                                "The given coverage is not an instance of 'TimeCoverage' or 'TimeLevelCoverage'")

    def write_variable_wind_stress(self):

        if (isinstance(self.coverage, TimeCoverage) or isinstance(self.coverage, TimeLevelCoverage)):

            if self.coverage.rank == 0:
                logging.info('[DefaultWriter] Gathering data from all processors')
                global_data = np.empty(
                    [self.coverage.get_t_size(type="target_global"),
                     2,
                     self.coverage.get_y_size(type="target_global"),
                     self.coverage.get_x_size(type="target_global")])
                global_data[:] = np.nan

            for time_index in range(0, self.coverage.get_t_size()):

                local_data = self.coverage.read_variable_wind_stress_at_time(time_index)

                if self.coverage.rank != 0:
                    self.coverage.comm.Send(np.ascontiguousarray(local_data), dest=0)
                else:
                    # Pour le proc n°1
                    global_data[self.coverage.map_mpi[self.coverage.rank]["dst_global_t"].start + time_index,
                                0,
                                self.coverage.map_mpi[self.coverage.rank]["dst_global_y"],
                                self.coverage.map_mpi[self.coverage.rank]["dst_global_x"]] = local_data[0]

                    global_data[self.coverage.map_mpi[self.coverage.rank]["dst_global_t"].start + time_index,
                                1,
                                self.coverage.map_mpi[self.coverage.rank]["dst_global_y"],
                                self.coverage.map_mpi[self.coverage.rank]["dst_global_x"]] = local_data[1]

                    # Pour les autres
                    for source in range(1, self.coverage.size):
                        recvbuf = np.empty([2, self.coverage.map_mpi[source]["dst_local_y_size"],
                                            self.coverage.map_mpi[source]["dst_local_x_size"]])
                        self.coverage.comm.Recv(recvbuf, source=source)

                        global_data[self.coverage.map_mpi[source]["dst_global_t"].start + time_index,
                                    0,
                                    self.coverage.map_mpi[source]["dst_global_y"],
                                    self.coverage.map_mpi[source]["dst_global_x"]] = recvbuf[0]

                        global_data[self.coverage.map_mpi[source]["dst_global_t"].start + time_index,
                                    1,
                                    self.coverage.map_mpi[source]["dst_global_y"],
                                    self.coverage.map_mpi[source]["dst_global_x"]] = recvbuf[1]

            self.coverage.comm.barrier()

            if self.coverage.rank == 0:

                logging.info('[DefaultWriter] Writing variable \'Wind Stress\'\'')

                for time_index in range(0, self.coverage.get_t_size(type="target_global")):
                    time = self.coverage.read_axis_t(type="target_global")[time_index]

                    logging.info('[DefaultWriter] Writing variable \'Wind Stress\' at time \'' + str(time) + '\'')

                    for vector in range(0, 2):
                        if vector == 0:

                            file = self.driver.Create(os.path.join(self.filename, time.strftime("%Y%m%d_%H%M%S") + "_" +
                                                                   VariableDefinition.VARIABLE_NAME[
                                                                       'eastward_wind_stress'] + ".tiff"),
                                                      int(self.rows), int(self.cols), 1, gdal.GDT_Float64)

                            # CRS info
                            proj = osr.SpatialReference()
                            proj.SetWellKnownGeogCS("EPSG:4326")
                            file.SetProjection(proj.ExportToWkt())
                            file.SetGeoTransform(self.geotransform)
                            file.GetRasterBand(1).WriteArray(global_data[time_index][0])

                        else:

                            file = self.driver.Create(os.path.join(self.filename, time.strftime("%Y%m%d_%H%M%S") + "_" +
                                                                   VariableDefinition.VARIABLE_NAME[
                                                                       'northward_wind_stress'] + ".tiff"),
                                                      int(self.rows), int(self.cols), 1, gdal.GDT_Float64)

                            # CRS info
                            proj = osr.SpatialReference()
                            proj.SetWellKnownGeogCS("EPSG:4326")
                            file.SetProjection(proj.ExportToWkt())
                            file.SetGeoTransform(self.geotransform)
                            file.GetRasterBand(1).WriteArray(global_data[time_index][1])

        else:
            raise CoverageError("DefaultWriter","The given coverage is not an instance of 'TimeCoverage' or 'TimeLevelCoverage'")

    #################
    # METEO
    # At 10 m
    #################
    def write_variable_wind_10m(self):

        if (isinstance(self.coverage, TimeCoverage) or isinstance(self.coverage, TimeLevelCoverage)):

            if self.coverage.rank == 0:
                logging.info('[DefaultWriter] Gathering data from all processors')
                global_data = np.empty(
                    [self.coverage.get_t_size(type="target_global"),
                     2,
                     self.coverage.get_y_size(type="target_global"),
                     self.coverage.get_x_size(type="target_global")])
                global_data[:] = np.nan

            for time_index in range(0, self.coverage.get_t_size()):

                local_data = self.coverage.read_variable_wind_10m_at_time(time_index)

                if self.coverage.rank != 0:
                    self.coverage.comm.Send(np.ascontiguousarray(local_data), dest=0)
                else:
                    # Pour le proc n°1
                    global_data[self.coverage.map_mpi[self.coverage.rank]["dst_global_t"].start + time_index,
                                0,
                                self.coverage.map_mpi[self.coverage.rank]["dst_global_y"],
                                self.coverage.map_mpi[self.coverage.rank]["dst_global_x"]] = local_data[0]

                    global_data[self.coverage.map_mpi[self.coverage.rank]["dst_global_t"].start + time_index,
                                1,
                                self.coverage.map_mpi[self.coverage.rank]["dst_global_y"],
                                self.coverage.map_mpi[self.coverage.rank]["dst_global_x"]] = local_data[1]

                    # Pour les autres
                    for source in range(1, self.coverage.size):
                        recvbuf = np.empty([2, self.coverage.map_mpi[source]["dst_local_y_size"],
                                            self.coverage.map_mpi[source]["dst_local_x_size"]])
                        self.coverage.comm.Recv(recvbuf, source=source)

                        global_data[self.coverage.map_mpi[source]["dst_global_t"].start + time_index,
                                    0,
                                    self.coverage.map_mpi[source]["dst_global_y"],
                                    self.coverage.map_mpi[source]["dst_global_x"]] = recvbuf[0]

                        global_data[self.coverage.map_mpi[source]["dst_global_t"].start + time_index,
                                    1,
                                    self.coverage.map_mpi[source]["dst_global_y"],
                                    self.coverage.map_mpi[source]["dst_global_x"]] = recvbuf[1]

            self.coverage.comm.barrier()

            if self.coverage.rank == 0:

                logging.info('[DefaultWriter] Writing variable \'Wind 10m\'')

                for time_index in range(0, self.coverage.get_t_size(type="target_global")):
                    time = self.coverage.read_axis_t(type="target_global")[time_index]

                    logging.info('[DefaultWriter] Writing variable \'Wind 10m\' at time \''+str(time)+'\'')

                    for vector in range(0, 2):
                        if vector == 0:

                            file = self.driver.Create(os.path.join(self.filename, time.strftime("%Y%m%d_%H%M%S") + "_" +
                                                                   VariableDefinition.VARIABLE_NAME[
                                                                       'eastward_wind_10m'] + ".tiff"),
                                                      int(self.rows), int(self.cols), 1, gdal.GDT_Float64)

                            # CRS info
                            proj = osr.SpatialReference()
                            proj.SetWellKnownGeogCS("EPSG:4326")
                            file.SetProjection(proj.ExportToWkt())
                            file.SetGeoTransform(self.geotransform)
                            file.GetRasterBand(1).WriteArray(global_data[time_index][0])

                        else:

                            file = self.driver.Create(os.path.join(self.filename, time.strftime("%Y%m%d_%H%M%S") + "_" +
                                                                   VariableDefinition.VARIABLE_NAME[
                                                                       'northward_wind_10m'] + ".tiff"),
                                                      int(self.rows), int(self.cols), 1, gdal.GDT_Float64)

                            # CRS info
                            proj = osr.SpatialReference()
                            proj.SetWellKnownGeogCS("EPSG:4326")
                            file.SetProjection(proj.ExportToWkt())
                            file.SetGeoTransform(self.geotransform)
                            file.GetRasterBand(1).WriteArray(global_data[time_index][1])
        else:
            raise CoverageError("DefaultWriter","The given coverage is not an instance of 'TimeCoverage' or 'TimeLevelCoverage'")

    def write_variable_wind_speed_10m(self):

        if (isinstance(self.coverage, TimeCoverage) or isinstance(self.coverage, TimeLevelCoverage)):

            if self.coverage.rank == 0:
                logging.info('[DefaultWriter] Gathering data from all processors')
                global_data = np.empty(
                    [self.coverage.get_t_size(type="target_global"),
                     self.coverage.get_y_size(type="target_global"),
                     self.coverage.get_x_size(type="target_global")])
                global_data[:] = np.nan

            for time_index in range(0, self.coverage.get_t_size()):

                local_data = self.coverage.read_variable_wind_speed_10m_at_time(
                    time_index)

                if self.coverage.rank != 0:
                    self.coverage.comm.Send(np.ascontiguousarray(local_data), dest=0)
                else:
                    # Pour le proc n°1
                    global_data[self.coverage.map_mpi[self.coverage.rank]["dst_global_t"].start + time_index,
                                self.coverage.map_mpi[self.coverage.rank]["dst_global_y"],
                                self.coverage.map_mpi[self.coverage.rank]["dst_global_x"]] = local_data

                    # Pour les autres
                    for source in range(1, self.coverage.size):
                        recvbuf = np.empty([self.coverage.map_mpi[source]["dst_local_y_size"],
                                            self.coverage.map_mpi[source]["dst_local_x_size"]])
                        self.coverage.comm.Recv(recvbuf, source=source)

                        global_data[self.coverage.map_mpi[source]["dst_global_t"].start + time_index,
                                    self.coverage.map_mpi[source]["dst_global_y"],
                                    self.coverage.map_mpi[source]["dst_global_x"]] = recvbuf

            self.coverage.comm.barrier()

            if self.coverage.rank == 0:

                logging.info('[DefaultWriter] Writing variable \'' + str(
                    VariableDefinition.LONG_NAME['wind_speed_10m']) + '\'')

                for time_index in range(0, self.coverage.get_t_size(type="target_global")):
                    time = self.coverage.read_axis_t(type="target_global")[time_index]

                    logging.debug('[DefaultWriter] Writing variable \'' + str(VariableDefinition.LONG_NAME[
                                                                                  'wind_speed_10m']) + '\' at time \'' + str(
                        time) + '\'')

                    file = self.driver.Create(os.path.join(self.filename, time.strftime("%Y%m%d_%H%M%S") + "_" +
                                                           VariableDefinition.VARIABLE_NAME[
                                                               'wind_speed_10m'] + ".tiff"),
                                              int(self.rows), int(self.cols), 1, gdal.GDT_Float64)

                    # CRS info
                    proj = osr.SpatialReference()
                    proj.SetWellKnownGeogCS("EPSG:4326")
                    file.SetProjection(proj.ExportToWkt())
                    file.SetGeoTransform(self.geotransform)
                    file.GetRasterBand(1).WriteArray(global_data[time_index])
        else:
            raise CoverageError("DefaultWriter","The given coverage is not an instance of 'TimeCoverage' or 'TimeLevelCoverage'")

    def write_variable_wind_to_direction_10m(self):

        if (isinstance(self.coverage, TimeCoverage) or isinstance(self.coverage, TimeLevelCoverage)):

            if self.coverage.rank == 0:
                logging.info('[DefaultWriter] Gathering data from all processors')
                global_data = np.empty(
                    [self.coverage.get_t_size(type="target_global"),
                     self.coverage.get_y_size(type="target_global"),
                     self.coverage.get_x_size(type="target_global")])
                global_data[:] = np.nan

            for time_index in range(0, self.coverage.get_t_size()):

                local_data = self.coverage.read_variable_wind_to_direction_10m_at_time(
                    time_index)

                if self.coverage.rank != 0:
                    self.coverage.comm.Send(np.ascontiguousarray(local_data), dest=0)
                else:
                    # Pour le proc n°1
                    global_data[self.coverage.map_mpi[self.coverage.rank]["dst_global_t"].start + time_index,
                                self.coverage.map_mpi[self.coverage.rank]["dst_global_y"],
                                self.coverage.map_mpi[self.coverage.rank]["dst_global_x"]] = local_data

                    # Pour les autres
                    for source in range(1, self.coverage.size):
                        recvbuf = np.empty([self.coverage.map_mpi[source]["dst_local_y_size"],
                                            self.coverage.map_mpi[source]["dst_local_x_size"]])
                        self.coverage.comm.Recv(recvbuf, source=source)

                        global_data[self.coverage.map_mpi[source]["dst_global_t"].start + time_index,
                                    self.coverage.map_mpi[source]["dst_global_y"],
                                    self.coverage.map_mpi[source]["dst_global_x"]] = recvbuf

            self.coverage.comm.barrier()

            if self.coverage.rank == 0:

                logging.info('[DefaultWriter] Writing variable \'' + str(
                    VariableDefinition.LONG_NAME['wind_to_direction_10m']) + '\'')

                for time_index in range(0, self.coverage.get_t_size(type="target_global")):
                    time = self.coverage.read_axis_t(type="target_global")[time_index]

                    logging.debug('[DefaultWriter] Writing variable \'' + str(VariableDefinition.LONG_NAME[
                                                                                  'wind_to_direction_10m']) + '\' at time \'' + str(
                        time) + '\'')

                    file = self.driver.Create(os.path.join(self.filename, time.strftime("%Y%m%d_%H%M%S") + "_" +
                                                           VariableDefinition.VARIABLE_NAME[
                                                               'wind_to_direction_10m'] + ".tiff"),
                                              int(self.rows), int(self.cols), 1, gdal.GDT_Float64)

                    # CRS info
                    proj = osr.SpatialReference()
                    proj.SetWellKnownGeogCS("EPSG:4326")
                    file.SetProjection(proj.ExportToWkt())
                    file.SetGeoTransform(self.geotransform)
                    file.GetRasterBand(1).WriteArray(global_data[time_index])
        else:
            raise CoverageError("DefaultWriter","The given coverage is not an instance of 'TimeCoverage' or 'TimeLevelCoverage'")

    def write_variable_wind_from_direction_10m(self):

        if (isinstance(self.coverage, TimeCoverage) or isinstance(self.coverage, TimeLevelCoverage)):

            if self.coverage.rank == 0:
                logging.info('[DefaultWriter] Gathering data from all processors')
                global_data = np.empty(
                    [self.coverage.get_t_size(type="target_global"),
                     self.coverage.get_y_size(type="target_global"),
                     self.coverage.get_x_size(type="target_global")])
                global_data[:] = np.nan

            for time_index in range(0, self.coverage.get_t_size()):

                local_data = self.coverage.read_variable_wind_from_direction_10m_at_time(
                    time_index)

                if self.coverage.rank != 0:
                    self.coverage.comm.Send(np.ascontiguousarray(local_data), dest=0)
                else:
                    # Pour le proc n°1
                    global_data[self.coverage.map_mpi[self.coverage.rank]["dst_global_t"].start + time_index,
                                self.coverage.map_mpi[self.coverage.rank]["dst_global_y"],
                                self.coverage.map_mpi[self.coverage.rank]["dst_global_x"]] = local_data

                    # Pour les autres
                    for source in range(1, self.coverage.size):
                        recvbuf = np.empty([self.coverage.map_mpi[source]["dst_local_y_size"],
                                            self.coverage.map_mpi[source]["dst_local_x_size"]])
                        self.coverage.comm.Recv(recvbuf, source=source)

                        global_data[self.coverage.map_mpi[source]["dst_global_t"].start + time_index,
                                    self.coverage.map_mpi[source]["dst_global_y"],
                                    self.coverage.map_mpi[source]["dst_global_x"]] = recvbuf

            self.coverage.comm.barrier()

            if self.coverage.rank == 0:

                logging.info('[DefaultWriter] Writing variable \'' + str(
                    VariableDefinition.LONG_NAME['wind_from_direction_10m']) + '\'')

                for time_index in range(0, self.coverage.get_t_size(type="target_global")):
                    time = self.coverage.read_axis_t(type="target_global")[time_index]

                    logging.debug('[DefaultWriter] Writing variable \'' + str(VariableDefinition.LONG_NAME[
                                                                                  'wind_from_direction_10m']) + '\' at time \'' + str(
                        time) + '\'')

                    file = self.driver.Create(os.path.join(self.filename, time.strftime("%Y%m%d_%H%M%S") + "_" +
                                                           VariableDefinition.VARIABLE_NAME[
                                                               'wind_from_direction_10m'] + ".tiff"),
                                              int(self.rows), int(self.cols), 1, gdal.GDT_Float64)

                    # CRS info
                    proj = osr.SpatialReference()
                    proj.SetWellKnownGeogCS("EPSG:4326")
                    file.SetProjection(proj.ExportToWkt())
                    file.SetGeoTransform(self.geotransform)
                    file.GetRasterBand(1).WriteArray(global_data[time_index])
        else:
            raise CoverageError("DefaultWriter","The given coverage is not an instance of 'TimeCoverage' or 'TimeLevelCoverage'")





