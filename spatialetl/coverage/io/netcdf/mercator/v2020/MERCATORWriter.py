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

import numpy as np
from mpi4py import MPI
from netCDF4 import Dataset
from netCDF4 import date2num
from numpy import int16, float32, float64

from spatialetl.coverage.LevelCoverage import LevelCoverage
from spatialetl.coverage.TimeCoverage import TimeCoverage
from spatialetl.coverage.TimeLevelCoverage import TimeLevelCoverage
from spatialetl.coverage.io.CoverageWriter import CoverageWriter
from spatialetl.exception.CoverageError import CoverageError
from spatialetl.utils.VariableDefinition import VariableDefinition
from spatialetl.utils.logger import logging


class DefaultWriter (CoverageWriter):

    def __init__(self,cov,myFile,mode='w'):
        CoverageWriter.__init__(self,cov,myFile);
        self.mode=mode
        format = 'NETCDF4_CLASSIC'

        if self.mode=='w':
            self.ncfile = Dataset(self.filename, 'w', parallel=True, comm=self.coverage.comm, info=MPI.Info(), format=format)
            self.ncfile.description = 'Generated with pySpatialETL'

            # dimensions
            self.ncfile.createDimension(VariableDefinition.VARIABLE_NAME['latitude'], self.coverage.get_y_size(type="target_global"))
            self.ncfile.createDimension(VariableDefinition.VARIABLE_NAME['longitude'], self.coverage.get_x_size(type="target_global"))

            if self.coverage.is_regular_grid()==True:

                # variables
                latitudes = self.ncfile.createVariable(VariableDefinition.VARIABLE_NAME['latitude'], float32, (VariableDefinition.VARIABLE_NAME['latitude'],))
                latitudes.long_name = VariableDefinition.LONG_NAME['latitude'] ;
                latitudes.standard_name = VariableDefinition.STANDARD_NAME['latitude'] ;
                latitudes.valid_min = -90.;
                latitudes.valid_max = 90. ;
                latitudes.axis = "Y" ;
                latitudes.units = VariableDefinition.CANONICAL_UNITS['latitude'];

                longitudes = self.ncfile.createVariable(VariableDefinition.VARIABLE_NAME['longitude'], float32, (VariableDefinition.VARIABLE_NAME['longitude'],))
                longitudes.long_name = VariableDefinition.LONG_NAME['longitude'] ;
                longitudes.standard_name = VariableDefinition.STANDARD_NAME['longitude'] ;
                longitudes.valid_min = -180. ;
                longitudes.valid_max = 180. ;
                longitudes.axis = "X" ;
                longitudes.units = VariableDefinition.CANONICAL_UNITS['longitude'];

                latitudes[self.coverage.map_mpi[self.coverage.rank]["dst_global_y"]] = self.coverage.read_axis_y()
                longitudes[self.coverage.map_mpi[self.coverage.rank]["dst_global_x"]] = self.coverage.read_axis_x()

                if(isinstance(self.coverage, TimeCoverage) or isinstance(self.coverage, TimeLevelCoverage)):

                    self.ncfile.createDimension(VariableDefinition.VARIABLE_NAME['time'], self.coverage.get_t_size(type="target_global"))
                    times = self.ncfile.createVariable(VariableDefinition.VARIABLE_NAME['time'], float64, (VariableDefinition.VARIABLE_NAME['time'],))
                    times.units= 'seconds since 1970-01-01 00:00:00'
                    times.calendar= 'gregorian'
                    times.standard_name= 'time'
                    times.axis='T'
                    times.conventions = "UTC time"

                    times[
                    self.coverage.map_mpi[self.coverage.rank]["dst_global_t"]] = date2num(self.coverage.read_axis_t(), units = times.units, calendar = times.calendar)

                if(isinstance(self.coverage, LevelCoverage) or isinstance(self.coverage, TimeLevelCoverage)):

                    if self.coverage.is_sigma_coordinate():
                        raise ValueError("This writer supports only Coverage with a regular vertical axis. Use 'resolution_z' on the Coverage to interpolate")

                    self.ncfile.createDimension(VariableDefinition.VARIABLE_NAME['depth'], self.coverage.get_z_size(type="target"))
                    levels = self.ncfile.createVariable(VariableDefinition.VARIABLE_NAME['depth'], float64, (VariableDefinition.VARIABLE_NAME['depth'],))
                    levels.standard_name= 'depth'
                    levels.long_name="Positive depth"
                    levels.axis='Z'
                    levels.units = VariableDefinition.CANONICAL_UNITS['depth'];

                    levels[:] = self.coverage.read_axis_z()
            else:
                raise ValueError("This writer supports only regular horizontal axis. Use 'resolution_x' and 'resolution_y' on the Coverage to interpolate")

        if self.mode=="a":

            self.ncfile = Dataset(self.filename, 'a', parallel=True, comm=self.coverage.comm, info=MPI.Info(),
                                  format=format)
            #self.ncfile.description = 'Generated with pySpatialETL'

            # dimensions
            if self.ncfile.dimensions['latitude'].size != self.coverage.get_y_size(type="target_global"):
                raise ValueError("Latitude dimension hasn't the same size than the Coverage. Unable to append the file.")

            if self.ncfile.dimensions['longitude'].size != self.coverage.get_x_size(type="target_global"):
                raise ValueError("Longitude dimension hasn't the same size than the Coverage. Unable to append the file.")

            if (isinstance(self.coverage, TimeCoverage) or isinstance(self.coverage, TimeLevelCoverage)):

                if self.ncfile.dimensions['time'].size != self.coverage.get_t_size(type="target_global"):
                    raise ValueError("Time dimension hasn't the same size than the Coverage. Unable to append the file.")

            if (isinstance(self.coverage, LevelCoverage) or isinstance(self.coverage, TimeLevelCoverage)):

                if self.coverage.is_sigma_coordinate():
                    raise ValueError("This writer supports only Coverage with a regular vertical axis.")

                if self.ncfile.dimensions['depth'].size != self.coverage.get_z_size(type="target"):
                    raise ValueError("Depth dimensions hasn't the same size than the Coverage. Unable to append the file.")

    def close(self):
        self.ncfile.close()

    # Variables
    def write_variable_2D_sea_binary_mask(self):

        if "zos" in self.ncfile.variables:
            var = self.ncfile.variables["zos"]
        else:
            var = self.ncfile.createVariable("zos", int16, (VariableDefinition.VARIABLE_NAME['latitude'],VariableDefinition.VARIABLE_NAME['longitude'],), fill_value=-9999)
        var.long_name = VariableDefinition.LONG_NAME['2d_sea_binary_mask']
        var.standard_name = VariableDefinition.STANDARD_NAME['2d_sea_binary_mask']
        var.units = VariableDefinition.CANONICAL_UNITS['2d_sea_binary_mask']
        var.comment = "1 = sea, 0 = land"

        if self.coverage.rank == 0:
            logging.info('[DefaultWriter] Writing variable \'' + str(VariableDefinition.LONG_NAME['2d_sea_binary_mask']) + '\'')

        var[
            self.coverage.map_mpi[self.coverage.rank]["dst_global_y"],
            self.coverage.map_mpi[self.coverage.rank]["dst_global_x"]
        ] = self.coverage.read_variable_2D_sea_binary_mask()

    #################
    # HYDRO
    # Sea Surface
    #################
    def write_variable_sea_surface_height_above_geoid(self):

        if (isinstance(self.coverage, TimeCoverage) or isinstance(self.coverage, TimeLevelCoverage)):

            if "zos" in self.ncfile.variables:
                var = self.ncfile.variables["zos"]
            else:
                var = self.ncfile.createVariable("zos", float32,
                                             (VariableDefinition.VARIABLE_NAME['time'], VariableDefinition.VARIABLE_NAME['latitude'], VariableDefinition.VARIABLE_NAME['longitude'],), fill_value=9.96921e+36)
            var.long_name = VariableDefinition.LONG_NAME['sea_surface_height_above_geoid']
            var.standard_name = VariableDefinition.STANDARD_NAME['sea_surface_height_above_geoid']
            var.units = VariableDefinition.CANONICAL_UNITS['sea_surface_height_above_geoid']

            if self.coverage.rank == 0:
                logging.info('[DefaultWriter] Writing variable \'' + str(VariableDefinition.LONG_NAME['sea_surface_height_above_geoid']) + '\'')

            time_index = 0
            for time in self.coverage.read_axis_t():
                logging.info('[DefaultWriter] Writing variable \'' + str(VariableDefinition.LONG_NAME['sea_surface_height_above_geoid']) + '\' at time \'' + str(time) + '\'')

                var[
                self.coverage.map_mpi[self.coverage.rank]["dst_global_t"].start+time_index:self.coverage.map_mpi[self.coverage.rank]["dst_global_t"].start+time_index+1,
                self.coverage.map_mpi[self.coverage.rank]["dst_global_y"],
                self.coverage.map_mpi[self.coverage.rank]["dst_global_x"]
                ] = self.coverage.read_variable_sea_surface_height_above_geoid_at_time(time)

                time_index += 1
        else:
            raise CoverageError("DefaultWriter","The given coverage is not an instance of 'TimeCoverage' or 'TimeLevelCoverage'")

    #################
    # HYDRO
    # 3D
    #################

    def write_variable_sea_water_temperature(self):

        if (isinstance(self.coverage, TimeLevelCoverage)):

            if "thetao" in self.ncfile.variables:
                var = self.ncfile.variables["thetao"]
            else:
                var = self.ncfile.createVariable("thetao", float32,
                                             (VariableDefinition.VARIABLE_NAME['time'], VariableDefinition.VARIABLE_NAME['depth'], VariableDefinition.VARIABLE_NAME['latitude'], VariableDefinition.VARIABLE_NAME['longitude'],), fill_value=9.96921e+36)
            var.long_name = VariableDefinition.LONG_NAME['sea_water_temperature']
            var.standard_name = VariableDefinition.STANDARD_NAME['sea_water_temperature']
            var.units = VariableDefinition.CANONICAL_UNITS['sea_water_temperature']

            if self.coverage.rank == 0:
                logging.info('[DefaultWriter] Writing variable \'' + str(VariableDefinition.LONG_NAME['sea_water_temperature']) + '\'')

            time_index = 0
            for time in self.coverage.read_axis_t():

                logging.info(
                    '[DefaultWriter] Writing variable \'' + str(VariableDefinition.LONG_NAME['sea_water_temperature']) + '\' at time \'' + str(time) + '\'')

                level_index = 0
                for level in self.coverage.read_axis_z():
                    # Pas d'interpolation temporelle donc on parcours les index du temps
                    var[
                    self.coverage.map_mpi[self.coverage.rank]["dst_global_t"].start + time_index:self.coverage.map_mpi[self.coverage.rank]["dst_global_t"].start + time_index + 1,
                    level_index:level_index + 1,
                    self.coverage.map_mpi[self.coverage.rank]["dst_global_y"],
                    self.coverage.map_mpi[self.coverage.rank]["dst_global_x"]
                    ] = self.coverage.read_variable_sea_water_temperature_at_time_and_depth(time, level)

                    level_index += 1

                time_index += 1
        else:
            raise CoverageError("DefaultWriter","The given coverage is not an instance of 'TimeLevelCoverage'")

    def write_variable_sea_water_salinity(self):

        if (isinstance(self.coverage, TimeLevelCoverage)):

            if "so" in self.ncfile.variables:
                var = self.ncfile.variables["so"]
            else:
                var = self.ncfile.createVariable("so", float32,
                                             (VariableDefinition.VARIABLE_NAME['time'], VariableDefinition.VARIABLE_NAME['depth'], VariableDefinition.VARIABLE_NAME['latitude'], VariableDefinition.VARIABLE_NAME['longitude'],), fill_value=9.96921e+36)
            var.long_name = VariableDefinition.LONG_NAME['sea_water_salinity']
            var.standard_name = VariableDefinition.STANDARD_NAME['sea_water_salinity']
            var.units = VariableDefinition.CANONICAL_UNITS['sea_water_salinity']

            if self.coverage.rank == 0:
                logging.info('[DefaultWriter] Writing variable \'' + str(VariableDefinition.LONG_NAME['sea_water_salinity']) + '\'')

            time_index = 0
            for time in self.coverage.read_axis_t():

                logging.info(
                '[DefaultWriter] Writing variable \'' + str(VariableDefinition.LONG_NAME['sea_water_salinity']) + '\' at time \'' + str(time) + '\'')

                level_index = 0
                for level in self.coverage.read_axis_z():
                    # Pas d'interpolation temporelle donc on parcours les index du temps
                    var[
                    self.coverage.map_mpi[self.coverage.rank]["dst_global_t"].start + time_index:self.coverage.map_mpi[self.coverage.rank]["dst_global_t"].start + time_index + 1,
                    level_index:level_index + 1,
                    self.coverage.map_mpi[self.coverage.rank]["dst_global_y"],
                    self.coverage.map_mpi[self.coverage.rank]["dst_global_x"]
                    ] = self.coverage.read_variable_sea_water_salinity_at_time_and_depth(time, level)

                    level_index += 1

                time_index += 1
        else:
            raise CoverageError("DefaultWriter","The given coverage is not an instance of 'TimeLevelCoverage'")

    def write_variable_baroclinic_sea_water_velocity(self):

        if (isinstance(self.coverage, TimeLevelCoverage)):

            if "uo" in self.ncfile.variables:
                ucomp = self.ncfile.variables["uo"]
            else:
                ucomp = self.ncfile.createVariable("uo",
                                              float32,
                                              (VariableDefinition.VARIABLE_NAME['time'],
                                               VariableDefinition.VARIABLE_NAME['depth'],
                                               VariableDefinition.VARIABLE_NAME['latitude'],
                                               VariableDefinition.VARIABLE_NAME['longitude'],), fill_value=9.96921e+36)
            ucomp.long_name = VariableDefinition.LONG_NAME['baroclinic_eastward_sea_water_velocity']
            ucomp.standard_name = VariableDefinition.STANDARD_NAME['baroclinic_eastward_sea_water_velocity']
            ucomp.units = VariableDefinition.CANONICAL_UNITS['baroclinic_eastward_sea_water_velocity']
            ucomp.comment = "cur=sqrt(U**2+V**2)";

            if "vo" in self.ncfile.variables:
                vcomp = self.ncfile.variables["vo"]
            else:
                vcomp = self.ncfile.createVariable("vo",
                                              float32,
                                              (VariableDefinition.VARIABLE_NAME['time'],
                                               VariableDefinition.VARIABLE_NAME['depth'],
                                               VariableDefinition.VARIABLE_NAME['latitude'],
                                               VariableDefinition.VARIABLE_NAME['longitude'],), fill_value=9.96921e+36)
            vcomp.long_name = VariableDefinition.LONG_NAME['baroclinic_northward_sea_water_velocity']
            vcomp.standard_name = VariableDefinition.STANDARD_NAME['baroclinic_northward_sea_water_velocity']
            vcomp.units = VariableDefinition.CANONICAL_UNITS['baroclinic_northward_sea_water_velocity']
            vcomp.comment = "cur=sqrt(U**2+V**2)";

            if self.coverage.rank == 0:
                logging.info('[DefaultWriter] Writing variable \'Baroclinic Sea Water Velocity\'\'')

            time_index = 0
            for time in self.coverage.read_axis_t():

                logging.info(
                    '[DefaultWriter] Writing variable \'Baroclinic Sea Water Velocity\' at time \'' + str(time) + '\'')

                level_index = 0
                for level in self.coverage.read_axis_z():
                    data_u,data_v = self.coverage.read_variable_baroclinic_sea_water_velocity_at_time_and_depth(time, level)

                    ucomp[
                    self.coverage.map_mpi[self.coverage.rank]["dst_global_t"].start + time_index:self.coverage.map_mpi[self.coverage.rank]["dst_global_t"].start + time_index + 1,
                    level_index:level_index + 1,
                    self.coverage.map_mpi[self.coverage.rank]["dst_global_y"],
                    self.coverage.map_mpi[self.coverage.rank]["dst_global_x"]
                    ] = data_u

                    vcomp[
                    self.coverage.map_mpi[self.coverage.rank]["dst_global_t"].start + time_index:self.coverage.map_mpi[self.coverage.rank]["dst_global_t"].start + time_index + 1,
                    level_index:level_index + 1,
                    self.coverage.map_mpi[self.coverage.rank]["dst_global_y"],
                    self.coverage.map_mpi[self.coverage.rank]["dst_global_x"]
                    ] = data_v
                    level_index += 1

                time_index += 1
        else:
            raise CoverageError("DefaultWriter","The given coverage is not an instance of 'TimeLevelCoverage'")

