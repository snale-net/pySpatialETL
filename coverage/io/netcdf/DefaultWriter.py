#! /usr/bin/env python2.7
# -*- coding: utf-8 -*-
#
# CoverageProcessing is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# any later version.
#
# CoverageProcessing is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# Author : Fabien Rétif - fabien.retif@zoho.com
#
from __future__ import division, print_function, absolute_import
from coverage.io.CoverageWriter import CoverageWriter
from utils.VariableDefinition import VariableDefinition
from netCDF4 import Dataset
from netCDF4 import date2num
from numpy import int16,float32,float64
import numpy as np
import logging
from coverage.TimeCoverage import TimeCoverage
from coverage.LevelCoverage import LevelCoverage
from coverage.TimeLevelCoverage import TimeLevelCoverage
from exception.CoverageError import CoverageError
from mpi4py import MPI

class DefaultWriter (CoverageWriter):

    def __init__(self,cov,myFile):
        CoverageWriter.__init__(self,cov,myFile);
        format = 'NETCDF3_CLASSIC'
        #format = 'NETCDF4_CLASSIC'
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

        else:

            latitudes = self.ncfile.createVariable(VariableDefinition.VARIABLE_NAME['latitude'], float32, (VariableDefinition.VARIABLE_NAME['latitude'],VariableDefinition.VARIABLE_NAME['longitude'],))
            latitudes.long_name = VariableDefinition.LONG_NAME['latitude'];
            latitudes.standard_name = VariableDefinition.STANDARD_NAME['latitude'];
            latitudes.valid_min = -90.;
            latitudes.valid_max = 90.;
            latitudes.axis = "Y";
            latitudes.units = VariableDefinition.CANONICAL_UNITS['latitude'];

            longitudes = self.ncfile.createVariable(VariableDefinition.VARIABLE_NAME['longitude'], float32, (VariableDefinition.VARIABLE_NAME['latitude'],VariableDefinition.VARIABLE_NAME['longitude'],))
            longitudes.long_name = VariableDefinition.LONG_NAME['longitude'];
            longitudes.standard_name = VariableDefinition.STANDARD_NAME['longitude'];
            longitudes.valid_min = -180.;
            longitudes.valid_max = 180.;
            longitudes.axis = "X";
            longitudes.units = VariableDefinition.CANONICAL_UNITS['longitude'];

            # data
            latitudes[
            self.coverage.map_mpi[self.coverage.rank]["dst_global_y"],self.coverage.map_mpi[self.coverage.rank]["dst_global_x"]] = self.coverage.read_axis_y()

            longitudes[
            self.coverage.map_mpi[self.coverage.rank]["dst_global_y"],self.coverage.map_mpi[self.coverage.rank]["dst_global_x"]] = self.coverage.read_axis_x()

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
                raise ValueError("This writer supports only Coverage with a regular vertical axis.")

            self.ncfile.createDimension(VariableDefinition.VARIABLE_NAME['depth'], self.coverage.get_z_size(type="target"))
            levels = self.ncfile.createVariable(VariableDefinition.VARIABLE_NAME['depth'], float64, (VariableDefinition.VARIABLE_NAME['depth'],))
            levels.standard_name= 'depth'
            levels.long_name="Positive depth"
            levels.axis='Z'
            levels.units = VariableDefinition.CANONICAL_UNITS['depth'];

            levels[:] = self.coverage.read_axis_z()

    def close(self):
        self.ncfile.close()

    # Variables
    def write_variable_mesh_size(self):

        var = self.ncfile.createVariable(VariableDefinition.VARIABLE_NAME['mesh_size'], float32, (VariableDefinition.VARIABLE_NAME['latitude'], VariableDefinition.VARIABLE_NAME['longitude'],), fill_value=9.96921e+36)
        var.long_name = VariableDefinition.LONG_NAME['mesh_size']
        var.standard_name = VariableDefinition.STANDARD_NAME['mesh_size']
        var.units = VariableDefinition.CANONICAL_UNITS['mesh_size']

        if self.coverage.rank == 0:
            logging.info('[DefaultWriter] Writing variable \'' + str(VariableDefinition.LONG_NAME['mesh_size']) + '\'')

        var[
            self.coverage.map_mpi[self.coverage.rank]["dst_global_y"],
            self.coverage.map_mpi[self.coverage.rank]["dst_global_x"]
        ] = self.coverage.read_variable_mesh_size()

    def write_variable_2D_sea_binary_mask(self):

        var = self.ncfile.createVariable(VariableDefinition.VARIABLE_NAME['2d_sea_binary_mask'], int16, (VariableDefinition.VARIABLE_NAME['latitude'],VariableDefinition.VARIABLE_NAME['longitude'],), fill_value=-9999)
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

    def write_variable_wet_binary_mask(self):

        if (isinstance(self.coverage, TimeCoverage) or isinstance(self.coverage, TimeLevelCoverage)):

            var = self.ncfile.createVariable(VariableDefinition.VARIABLE_NAME['wet_binary_mask'], float32, (VariableDefinition.VARIABLE_NAME['time'], VariableDefinition.VARIABLE_NAME['latitude'], VariableDefinition.VARIABLE_NAME['longitude'],),
                                             fill_value=9.96921e+36)
            var.long_name = VariableDefinition.LONG_NAME['wet_binary_mask']
            var.standard_name = VariableDefinition.STANDARD_NAME['wet_binary_mask']
            var.units = VariableDefinition.CANONICAL_UNITS['wet_binary_mask']
            var.comment = "1 = sea, 0 = land"

            if self.coverage.rank == 0:
                logging.info('[DefaultWriter] Writing variable \'' + str(VariableDefinition.VARIABLE_NAME['wet_binary_mask']) + '\'')

            time_index = 0
            for time in self.coverage.read_axis_t():
                logging.debug(
                    '[DefaultWriter] Writing variable \'' + str(VariableDefinition.VARIABLE_NAME['wet_binary_mask']) + '\' at time \'' + str(time) + '\'')
                var[
                self.coverage.map_mpi[self.coverage.rank]["dst_global_t"].start+time_index:self.coverage.map_mpi[self.coverage.rank]["dst_global_t"].start+time_index+1,
                self.coverage.map_mpi[self.coverage.rank]["dst_global_y"],
                self.coverage.map_mpi[self.coverage.rank]["dst_global_x"]
                ] = self.coverage.read_variable_2D_wet_binary_mask_at_time(time)

                time_index += 1

        else:
            raise CoverageError("DefaultWriter","The given coverage is not an instance of 'TimeCoverage' or 'TimeLevelCoverage'")

    #################
    # HYDRO
    # 2D
    #################
    def write_variable_bathymetry(self):

        var = self.ncfile.createVariable(VariableDefinition.VARIABLE_NAME['bathymetry'], float32, (VariableDefinition.VARIABLE_NAME['latitude'], VariableDefinition.VARIABLE_NAME['longitude'],),fill_value=9.96921e+36)
        var.long_name = VariableDefinition.LONG_NAME['bathymetry']
        var.standard_name = VariableDefinition.STANDARD_NAME['bathymetry']
        var.units = VariableDefinition.CANONICAL_UNITS['bathymetry']

        if self.coverage.rank == 0:
            logging.info('[DefaultWriter] Writing variable \'' + str(VariableDefinition.LONG_NAME['bathymetry']) + '\'')

        var[
        self.coverage.map_mpi[self.coverage.rank]["dst_global_y"],
        self.coverage.map_mpi[self.coverage.rank]["dst_global_x"]
        ] = self.coverage.read_variable_bathymetry()

    def write_variable_barotropic_sea_water_velocity(self):

        if (isinstance(self.coverage, TimeCoverage) or isinstance(self.coverage, TimeLevelCoverage)):

            ucur = self.ncfile.createVariable(VariableDefinition.VARIABLE_NAME['barotropic_eastward_sea_water_velocity'],
                                              float32,
                                              (VariableDefinition.VARIABLE_NAME['time'],
                                               VariableDefinition.VARIABLE_NAME['latitude'],
                                               VariableDefinition.VARIABLE_NAME['longitude'],), fill_value=9.96921e+36)
            ucur.long_name = VariableDefinition.LONG_NAME['barotropic_eastward_sea_water_velocity']
            ucur.standard_name = VariableDefinition.STANDARD_NAME['barotropic_eastward_sea_water_velocity']
            ucur.units = VariableDefinition.CANONICAL_UNITS['barotropic_eastward_sea_water_velocity']
            ucur.comment = "cur=sqrt(U**2+V**2)";

            vcur = self.ncfile.createVariable(VariableDefinition.VARIABLE_NAME['barotropic_northward_sea_water_velocity'],
                                              float32,
                                              (VariableDefinition.VARIABLE_NAME['time'],
                                               VariableDefinition.VARIABLE_NAME['latitude'],
                                               VariableDefinition.VARIABLE_NAME['longitude'],), fill_value=9.96921e+36)
            vcur.long_name = VariableDefinition.LONG_NAME['barotropic_northward_sea_water_velocity']
            vcur.standard_name = VariableDefinition.STANDARD_NAME['barotropic_northward_sea_water_velocity']
            vcur.units = VariableDefinition.CANONICAL_UNITS['barotropic_northward_sea_water_velocity']
            vcur.comment = "cur=sqrt(U**2+V**2)";

            if self.coverage.rank == 0:
                logging.info('[DefaultWriter] Writing variable \'Barotropic Sea Water Velocity\'\'')

            time_index = 0
            for time in self.coverage.read_axis_t():
                logging.debug(
                    '[DefaultWriter] Writing variable \'Barotropic Sea Water Velocity\' at time \'' + str(time) + '\'')

                cur = self.coverage.read_variable_barotropic_sea_water_velocity_at_time(time_index)

                ucur[ self.coverage.map_mpi[self.coverage.rank]["dst_global_t"].start+time_index:self.coverage.map_mpi[self.coverage.rank]["dst_global_t"].start+time_index+1,
                self.coverage.map_mpi[self.coverage.rank]["dst_global_y"],
                self.coverage.map_mpi[self.coverage.rank]["dst_global_x"]] = cur[0]

                vcur[ self.coverage.map_mpi[self.coverage.rank]["dst_global_t"].start+time_index:self.coverage.map_mpi[self.coverage.rank]["dst_global_t"].start+time_index+1,
                self.coverage.map_mpi[self.coverage.rank]["dst_global_y"],
                self.coverage.map_mpi[self.coverage.rank]["dst_global_x"]] = cur[1]

                time_index += 1

        else:
            raise CoverageError("DefaultWriter","The given coverage is not an instance of 'TimeCoverage' or 'TimeLevelCoverage'")

    #################
    # HYDRO
    # Sea Surface
    #################
    def write_variable_sea_surface_height_above_mean_sea_level(self):

        if (isinstance(self.coverage, TimeCoverage) or isinstance(self.coverage, TimeLevelCoverage)):
            
            var = self.ncfile.createVariable(VariableDefinition.VARIABLE_NAME['sea_surface_height_above_mean_sea_level'], float32, (VariableDefinition.VARIABLE_NAME['time'], VariableDefinition.VARIABLE_NAME['latitude'], VariableDefinition.VARIABLE_NAME['longitude'],),fill_value=9.96921e+36)
            var.long_name = VariableDefinition.LONG_NAME['sea_surface_height_above_mean_sea_level']
            var.standard_name = VariableDefinition.STANDARD_NAME['sea_surface_height_above_mean_sea_level']
            var.units = VariableDefinition.CANONICAL_UNITS['sea_surface_height_above_mean_sea_level']

            if self.coverage.rank == 0:
                logging.info('[DefaultWriter] Writing variable \''+str(VariableDefinition.LONG_NAME['sea_surface_height_above_mean_sea_level'])+'\'')

            time_index=0
            for time in self.coverage.read_axis_t():

                logging.info('[DefaultWriter] Writing variable \''+str(VariableDefinition.LONG_NAME['sea_surface_height_above_mean_sea_level'])+'\' at time \''+str(time)+'\'')

                var[
                self.coverage.map_mpi[self.coverage.rank]["dst_global_t"].start+time_index:self.coverage.map_mpi[self.coverage.rank]["dst_global_t"].start+time_index+1,
                self.coverage.map_mpi[self.coverage.rank]["dst_global_y"],
                self.coverage.map_mpi[self.coverage.rank]["dst_global_x"]
                ] = self.coverage.read_variable_sea_surface_height_above_mean_sea_level_at_time(time)

                time_index += 1
        else:
            raise CoverageError("DefaultWriter","The given coverage is not an instance of 'TimeCoverage' or 'TimeLevelCoverage'")

    def write_variable_sea_surface_height_above_geoid(self):

        if (isinstance(self.coverage, TimeCoverage) or isinstance(self.coverage, TimeLevelCoverage)):

            var = self.ncfile.createVariable(VariableDefinition.VARIABLE_NAME['sea_surface_height_above_geoid'], float32,
                                             (VariableDefinition.VARIABLE_NAME['time'], VariableDefinition.VARIABLE_NAME['latitude'], VariableDefinition.VARIABLE_NAME['longitude'],), fill_value=9.96921e+36)
            var.long_name = VariableDefinition.LONG_NAME['sea_surface_height_above_geoid']
            var.standard_name = VariableDefinition.STANDARD_NAME['sea_surface_height_above_geoid']
            var.units = VariableDefinition.CANONICAL_UNITS['sea_surface_height_above_geoid']

            if self.coverage.rank == 0:
                logging.info('[DefaultWriter] Writing variable \'' + str(VariableDefinition.LONG_NAME['sea_surface_height_above_geoid']) + '\'')

            time_index = 0
            for time in self.coverage.read_axis_t():
                logging.debug('[DefaultWriter] Writing variable \'' + str(VariableDefinition.LONG_NAME['sea_surface_height_above_geoid']) + '\' at time \'' + str(time) + '\'')

                var[
                self.coverage.map_mpi[self.coverage.rank]["dst_global_t"].start+time_index:self.coverage.map_mpi[self.coverage.rank]["dst_global_t"].start+time_index+1,
                self.coverage.map_mpi[self.coverage.rank]["dst_global_y"],
                self.coverage.map_mpi[self.coverage.rank]["dst_global_x"]
                ] = self.coverage.read_variable_sea_surface_height_above_geoid_at_time(time_index)

                time_index += 1

        else:
            raise CoverageError("DefaultWriter","The given coverage is not an instance of 'TimeCoverage' or 'TimeLevelCoverage'")

    def write_variable_sea_surface_temperature(self):

        if (isinstance(self.coverage, TimeCoverage) or isinstance(self.coverage, TimeLevelCoverage)):

            var = self.ncfile.createVariable(VariableDefinition.VARIABLE_NAME['sea_surface_temperature'], float32, (
            VariableDefinition.VARIABLE_NAME['time'], VariableDefinition.VARIABLE_NAME['latitude'],
            VariableDefinition.VARIABLE_NAME['longitude'],),
                                             fill_value=9.96921e+36)
            var.long_name = VariableDefinition.LONG_NAME['sea_surface_temperature']
            var.standard_name = VariableDefinition.STANDARD_NAME['sea_surface_temperature']
            var.units = VariableDefinition.CANONICAL_UNITS['sea_surface_temperature']

            if self.coverage.rank == 0:
                logging.info('[DefaultWriter] Writing variable \'' + str(VariableDefinition.LONG_NAME['sea_surface_temperature']) + '\'')

            time_index = 0
            for time in self.coverage.read_axis_t():
                logging.debug(
                    '[DefaultWriter] Writing variable \'' + str(
                        VariableDefinition.LONG_NAME['sea_surface_temperature']) + '\' at time \'' + str(time) + '\'')

                var[
                self.coverage.map_mpi[self.coverage.rank]["dst_global_t"].start+time_index:self.coverage.map_mpi[self.coverage.rank]["dst_global_t"].start+time_index+1,
                self.coverage.map_mpi[self.coverage.rank]["dst_global_y"],
                self.coverage.map_mpi[self.coverage.rank]["dst_global_x"]
                ] = self.coverage.read_variable_sea_surface_temperature_at_time(time_index)

                time_index += 1

        else:
            raise CoverageError("DefaultWriter","The given coverage is not an instance of 'TimeCoverage' or 'TimeLevelCoverage'")

    def write_variable_sea_surface_salinity(self):

        if (isinstance(self.coverage, TimeCoverage) or isinstance(self.coverage, TimeLevelCoverage)):

            var = self.ncfile.createVariable(VariableDefinition.VARIABLE_NAME['sea_surface_salinity'], float32, (
            VariableDefinition.VARIABLE_NAME['time'], VariableDefinition.VARIABLE_NAME['latitude'],
            VariableDefinition.VARIABLE_NAME['longitude'],),
                                             fill_value=9.96921e+36)
            var.long_name = VariableDefinition.LONG_NAME['sea_surface_salinity']
            var.standard_name = VariableDefinition.STANDARD_NAME['sea_surface_salinity']
            var.units = VariableDefinition.CANONICAL_UNITS['sea_surface_salinity']

            if self.coverage.rank == 0:
                logging.info('[DefaultWriter] Writing variable \'' + str(VariableDefinition.LONG_NAME['sea_surface_salinity']) + '\'')

            time_index = 0
            for time in self.coverage.read_axis_t():
                logging.debug(
                    '[DefaultWriter] Writing variable \'' + str(
                        VariableDefinition.LONG_NAME['sea_surface_salinity']) + '\' at time \'' + str(time) + '\'')

                var[
                self.coverage.map_mpi[self.coverage.rank]["dst_global_t"].start+time_index:self.coverage.map_mpi[self.coverage.rank]["dst_global_t"].start+time_index+1,
                self.coverage.map_mpi[self.coverage.rank]["dst_global_y"],
                self.coverage.map_mpi[self.coverage.rank]["dst_global_x"]
                ] = self.coverage.read_variable_sea_surface_salinity_at_time(time_index)

                time_index += 1

        else:
            raise CoverageError("DefaultWriter","The given coverage is not an instance of 'TimeCoverage' or 'TimeLevelCoverage'")

    #################
    # HYDRO
    # Ground level
    #################
    def write_variable_sea_water_temperature_at_ground_level(self):

        if (isinstance(self.coverage, TimeCoverage) or isinstance(self.coverage, TimeLevelCoverage)):

            var = self.ncfile.createVariable(VariableDefinition.VARIABLE_NAME['sea_water_temperature_at_ground_level'], float32, (
            VariableDefinition.VARIABLE_NAME['time'], VariableDefinition.VARIABLE_NAME['latitude'],
            VariableDefinition.VARIABLE_NAME['longitude'],),
                                             fill_value=9.96921e+36)
            var.long_name = VariableDefinition.LONG_NAME['sea_water_temperature_at_ground_level']
            var.standard_name = VariableDefinition.STANDARD_NAME['sea_water_temperature_at_ground_level']
            var.units = VariableDefinition.CANONICAL_UNITS['sea_water_temperature_at_ground_level']

            if self.coverage.rank == 0:
                logging.info('[DefaultWriter] Writing variable \'' + str(VariableDefinition.LONG_NAME['sea_water_temperature_at_ground_level']) + '\'')

            time_index = 0
            for time in self.coverage.read_axis_t():
                logging.debug(
                    '[DefaultWriter] Writing variable \'' + str(
                        VariableDefinition.LONG_NAME['sea_water_temperature_at_ground_level']) + '\' at time \'' + str(time) + '\'')

                var[
                self.coverage.map_mpi[self.coverage.rank]["dst_global_t"].start+time_index:self.coverage.map_mpi[self.coverage.rank]["dst_global_t"].start+time_index+1,
                self.coverage.map_mpi[self.coverage.rank]["dst_global_y"],
                self.coverage.map_mpi[self.coverage.rank]["dst_global_x"]
                ] = self.coverage.read_variable_sea_water_temperature_at_ground_level_at_time(time_index)

                time_index += 1

        else:
            raise CoverageError("DefaultWriter","The given coverage is not an instance of 'TimeCoverage' or 'TimeLevelCoverage'")

    def write_variable_sea_water_salinity_at_ground_level(self):

        if (isinstance(self.coverage, TimeCoverage) or isinstance(self.coverage, TimeLevelCoverage)):

            var = self.ncfile.createVariable(VariableDefinition.VARIABLE_NAME['sea_water_salinity_at_ground_level'], float32, (
            VariableDefinition.VARIABLE_NAME['time'], VariableDefinition.VARIABLE_NAME['latitude'],
            VariableDefinition.VARIABLE_NAME['longitude'],),
                                             fill_value=9.96921e+36)
            var.long_name = VariableDefinition.LONG_NAME['sea_water_salinity_at_ground_level']
            var.standard_name = VariableDefinition.STANDARD_NAME['sea_water_salinity_at_ground_level']
            var.units = VariableDefinition.CANONICAL_UNITS['sea_water_salinity_at_ground_level']

            if self.coverage.rank == 0:
                logging.info('[DefaultWriter] Writing variable \'' + str(VariableDefinition.LONG_NAME['sea_water_salinity_at_ground_level']) + '\'')

            time_index = 0
            for time in self.coverage.read_axis_t():
                logging.debug(
                    '[DefaultWriter] Writing variable \'' + str(
                        VariableDefinition.LONG_NAME['sea_water_salinity_at_ground_level']) + '\' at time \'' + str(time) + '\'')

                var[
                self.coverage.map_mpi[self.coverage.rank]["dst_global_t"].start+time_index:self.coverage.map_mpi[self.coverage.rank]["dst_global_t"].start+time_index+1,
                self.coverage.map_mpi[self.coverage.rank]["dst_global_y"],
                self.coverage.map_mpi[self.coverage.rank]["dst_global_x"]
                ] = self.coverage.read_variable_sea_water_salinity_at_ground_level_at_time(time_index)

                time_index += 1

        else:
            raise CoverageError("DefaultWriter","The given coverage is not an instance of 'TimeCoverage' or 'TimeLevelCoverage'")

    #################
    # HYDRO
    # 3D
    #################
    def write_variable_sea_water_temperature(self):

        if (isinstance(self.coverage, TimeLevelCoverage)):

            var = self.ncfile.createVariable(VariableDefinition.VARIABLE_NAME['sea_water_temperature'], float32,
                                             (VariableDefinition.VARIABLE_NAME['time'], VariableDefinition.VARIABLE_NAME['depth'], VariableDefinition.VARIABLE_NAME['latitude'], VariableDefinition.VARIABLE_NAME['longitude'],), fill_value=9.96921e+36)
            var.long_name = VariableDefinition.LONG_NAME['sea_water_temperature']
            var.standard_name = VariableDefinition.STANDARD_NAME['sea_water_temperature']
            var.units = VariableDefinition.CANONICAL_UNITS['sea_water_temperature']

            if self.coverage.rank == 0:
                logging.info('[DefaultWriter] Writing variable \'' + str(VariableDefinition.LONG_NAME['sea_water_temperature']) + '\'')

            time_index = 0
            for time in self.coverage.read_axis_t():

                logging.debug(
                    '[DefaultWriter] Writing variable \'' + str(VariableDefinition.LONG_NAME['sea_water_temperature']) + '\' at time \'' + str(time) + '\'')

                level_index = 0
                for level in self.coverage.read_axis_z():
                    # Pas d'interpolation temporelle donc on parcours les index du temps
                    var[
                    self.coverage.map_mpi[self.coverage.rank]["dst_global_t"].start + time_index:self.coverage.map_mpi[self.coverage.rank]["dst_global_t"].start + time_index + 1,
                    level_index:level_index + 1,
                    self.coverage.map_mpi[self.coverage.rank]["dst_global_y"],
                    self.coverage.map_mpi[self.coverage.rank]["dst_global_x"]
                    ] = self.coverage.read_variable_sea_water_temperature_at_time_and_depth(time_index, level)

                    level_index += 1

                time_index += 1

        else:
            raise CoverageError("DefaultWriter","The given coverage is not an instance of 'TimeLevelCoverage'")

    def write_variable_sea_water_salinity(self):

        if (isinstance(self.coverage, TimeLevelCoverage)):

            var = self.ncfile.createVariable(VariableDefinition.VARIABLE_NAME['sea_water_salinity'], float32,
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
                    ] = self.coverage.read_variable_sea_water_salinity_at_time_and_depth(time_index, level)

                    level_index += 1

                time_index += 1

        else:
            raise CoverageError("DefaultWriter","The given coverage is not an instance of 'TimeLevelCoverage'")

    def write_variable_baroclinic_sea_water_velocity(self):

        if (isinstance(self.coverage, TimeLevelCoverage)):

            ucur = self.ncfile.createVariable(VariableDefinition.VARIABLE_NAME['baroclinic_eastward_sea_water_velocity'],
                                              float32,
                                              (VariableDefinition.VARIABLE_NAME['time'],
                                               VariableDefinition.VARIABLE_NAME['depth'],
                                               VariableDefinition.VARIABLE_NAME['latitude'],
                                               VariableDefinition.VARIABLE_NAME['longitude'],), fill_value=9.96921e+36)
            ucur.long_name = VariableDefinition.LONG_NAME['baroclinic_eastward_sea_water_velocity']
            ucur.standard_name = VariableDefinition.STANDARD_NAME['baroclinic_eastward_sea_water_velocity']
            ucur.units = VariableDefinition.CANONICAL_UNITS['baroclinic_eastward_sea_water_velocity']
            ucur.comment = "cur=sqrt(U**2+V**2)";

            vcur = self.ncfile.createVariable(VariableDefinition.VARIABLE_NAME['baroclinic_northward_sea_water_velocity'],
                                              float32,
                                              (VariableDefinition.VARIABLE_NAME['time'],
                                               VariableDefinition.VARIABLE_NAME['depth'],
                                               VariableDefinition.VARIABLE_NAME['latitude'],
                                               VariableDefinition.VARIABLE_NAME['longitude'],), fill_value=9.96921e+36)
            vcur.long_name = VariableDefinition.LONG_NAME['baroclinic_northward_sea_water_velocity']
            vcur.standard_name = VariableDefinition.STANDARD_NAME['baroclinic_northward_sea_water_velocity']
            vcur.units = VariableDefinition.CANONICAL_UNITS['baroclinic_northward_sea_water_velocity']
            vcur.comment = "cur=sqrt(U**2+V**2)";

            if self.coverage.rank == 0:
                logging.info('[DefaultWriter] Writing variable \'Baroclinic Sea Water Velocity\'\'')

            time_index = 0
            for time in self.coverage.read_axis_t():

                logging.debug(
                    '[DefaultWriter] Writing variable \'Baroclinic Sea Water Velocity\' at time \'' + str(time) + '\'')

                level_index = 0
                for level in self.coverage.read_axis_z():
                    # Pas d'interpolation temporelle donc on parcours les index du temps
                    cur = self.coverage.read_variable_baroclinic_sea_water_velocity_at_time_and_depth(time_index, level)

                    ucur[
                    self.coverage.map_mpi[self.coverage.rank]["dst_global_t"].start + time_index:self.coverage.map_mpi[self.coverage.rank]["dst_global_t"].start + time_index + 1,
                    level_index:level_index + 1,
                    self.coverage.map_mpi[self.coverage.rank]["dst_global_y"],
                    self.coverage.map_mpi[self.coverage.rank]["dst_global_x"]
                    ] = cur[0]

                    vcur[
                    self.coverage.map_mpi[self.coverage.rank]["dst_global_t"].start + time_index:self.coverage.map_mpi[self.coverage.rank]["dst_global_t"].start + time_index + 1,
                    level_index:level_index + 1,
                    self.coverage.map_mpi[self.coverage.rank]["dst_global_y"],
                    self.coverage.map_mpi[self.coverage.rank]["dst_global_x"]
                    ] = cur[1]
                    level_index += 1

                time_index += 1

        else:
            raise CoverageError("DefaultWriter","The given coverage is not an instance of 'TimeLevelCoverage'")

    #################
    # WAVES
    # Sea Surface
    #################
    def write_variable_sea_surface_wave_significant_height(self):

        if (isinstance(self.coverage, TimeCoverage) or isinstance(self.coverage, TimeLevelCoverage)):

            var = self.ncfile.createVariable(VariableDefinition.VARIABLE_NAME['sea_surface_wave_significant_height'], float32, (VariableDefinition.VARIABLE_NAME['time'], VariableDefinition.VARIABLE_NAME['latitude'], VariableDefinition.VARIABLE_NAME['longitude'],),fill_value=9.96921e+36)
            var.long_name = VariableDefinition.LONG_NAME['sea_surface_wave_significant_height']
            var.standard_name = VariableDefinition.STANDARD_NAME['sea_surface_wave_significant_height']
            var.units = VariableDefinition.CANONICAL_UNITS['sea_surface_wave_significant_height']

            if self.coverage.rank == 0:
                logging.info('[DefaultWriter] Writing variable \'' + str(VariableDefinition.LONG_NAME['sea_surface_wave_significant_height']) + '\'')

            time_index=0
            for time in self.coverage.read_axis_t():
                logging.debug(
                    '[DefaultWriter] Writing variable \'' + str(VariableDefinition.LONG_NAME['sea_surface_wave_significant_height']) + '\' at time \'' + str(time) + '\'')

                var[
                self.coverage.map_mpi[self.coverage.rank]["dst_global_t"].start+time_index:self.coverage.map_mpi[self.coverage.rank]["dst_global_t"].start+time_index+1,
                self.coverage.map_mpi[self.coverage.rank]["dst_global_y"],
                self.coverage.map_mpi[self.coverage.rank]["dst_global_x"]
                ] = self.coverage.read_variable_sea_surface_wave_significant_height_at_time(time_index)

                time_index += 1

        else:
            raise CoverageError("DefaultWriter","The given coverage is not an instance of 'TimeCoverage' or 'TimeLevelCoverage'")

    def write_variable_sea_surface_wave_breaking_height(self):

        if (isinstance(self.coverage, TimeCoverage) or isinstance(self.coverage, TimeLevelCoverage)):

            var = self.ncfile.createVariable(VariableDefinition.VARIABLE_NAME['sea_surface_wave_breaking_height'], float32,
                                             (VariableDefinition.VARIABLE_NAME['time'], VariableDefinition.VARIABLE_NAME['latitude'], VariableDefinition.VARIABLE_NAME['longitude'],), fill_value=9.96921e+36)
            var.long_name = VariableDefinition.LONG_NAME['sea_surface_wave_breaking_height']
            var.standard_name = VariableDefinition.STANDARD_NAME['sea_surface_wave_breaking_height']
            var.units =VariableDefinition.CANONICAL_UNITS['sea_surface_wave_breaking_height']

            if self.coverage.rank == 0:
                logging.info('[DefaultWriter] Writing variable \'' + str(VariableDefinition.LONG_NAME['sea_surface_wave_breaking_height']) + '\'')

            time_index = 0
            for time in self.coverage.read_axis_t():
                logging.debug(
                    '[DefaultWriter] Writing variable \'' + str(VariableDefinition.LONG_NAME['sea_surface_wave_breaking_height']) + '\' at time \'' + str(time) + '\'')

                var[
                self.coverage.map_mpi[self.coverage.rank]["dst_global_t"].start+time_index:self.coverage.map_mpi[self.coverage.rank]["dst_global_t"].start+time_index+1,
                self.coverage.map_mpi[self.coverage.rank]["dst_global_y"],
                self.coverage.map_mpi[self.coverage.rank]["dst_global_x"]
                ] = self.coverage.read_variable_sea_surface_wave_breaking_height_at_time(time_index)

                time_index += 1

        else:
            raise CoverageError("DefaultWriter","The given coverage is not an instance of 'TimeCoverage' or 'TimeLevelCoverage'")

    def write_variable_sea_surface_wave_mean_period(self):

        if (isinstance(self.coverage, TimeCoverage) or isinstance(self.coverage, TimeLevelCoverage)):

            var = self.ncfile.createVariable(VariableDefinition.VARIABLE_NAME['sea_surface_wave_mean_period'], float32,
                                             (VariableDefinition.VARIABLE_NAME['time'], VariableDefinition.VARIABLE_NAME['latitude'], VariableDefinition.VARIABLE_NAME['longitude'],), fill_value=9.96921e+36)
            var.long_name = VariableDefinition.LONG_NAME['sea_surface_wave_mean_period']
            var.standard_name = VariableDefinition.STANDARD_NAME['sea_surface_wave_mean_period']
            var.units = VariableDefinition.CANONICAL_UNITS['sea_surface_wave_mean_period']

            if self.coverage.rank == 0:
                logging.info('[DefaultWriter] Writing variable \'' + str(VariableDefinition.LONG_NAME['sea_surface_wave_mean_period']) + '\'')

            time_index = 0
            for time in self.coverage.read_axis_t():
                logging.debug(
                    '[DefaultWriter] Writing variable \'' + str(VariableDefinition.LONG_NAME['sea_surface_wave_mean_period']) + '\' at time \'' + str(time) + '\'')

                var[
                self.coverage.map_mpi[self.coverage.rank]["dst_global_t"].start+time_index:self.coverage.map_mpi[self.coverage.rank]["dst_global_t"].start+time_index+1,
                self.coverage.map_mpi[self.coverage.rank]["dst_global_y"],
                self.coverage.map_mpi[self.coverage.rank]["dst_global_x"]
                ] = self.coverage.read_variable_sea_surface_wave_mean_period_at_time(
                    time_index)

                time_index += 1

        else:
            raise CoverageError("DefaultWriter","The given coverage is not an instance of 'TimeCoverage' or 'TimeLevelCoverage'")

    def write_variable_sea_surface_wave_peak_period(self):

        if (isinstance(self.coverage, TimeCoverage) or isinstance(self.coverage, TimeLevelCoverage)):

            var = self.ncfile.createVariable(VariableDefinition.VARIABLE_NAME['sea_surface_wave_peak_period'], float32,
                                             (VariableDefinition.VARIABLE_NAME['time'], VariableDefinition.VARIABLE_NAME['latitude'], VariableDefinition.VARIABLE_NAME['longitude'],), fill_value=9.96921e+36)
            var.long_name = VariableDefinition.LONG_NAME['sea_surface_wave_peak_period']
            var.standard_name = VariableDefinition.STANDARD_NAME['sea_surface_wave_peak_period']
            var.units = VariableDefinition.CANONICAL_UNITS['sea_surface_wave_peak_period']

            if self.coverage.rank == 0:
                logging.info('[DefaultWriter] Writing variable \'' + str(VariableDefinition.LONG_NAME['sea_surface_wave_peak_period']) + '\'')

            time_index = 0
            for time in self.coverage.read_axis_t():
                logging.debug(
                    '[DefaultWriter] Writing variable \'' + str(VariableDefinition.LONG_NAME['sea_surface_wave_peak_period']) + '\' at time \'' + str(time) + '\'')

                var[
                self.coverage.map_mpi[self.coverage.rank]["dst_global_t"].start+time_index:self.coverage.map_mpi[self.coverage.rank]["dst_global_t"].start+time_index+1,
                self.coverage.map_mpi[self.coverage.rank]["dst_global_y"],
                self.coverage.map_mpi[self.coverage.rank]["dst_global_x"]
                ] = self.coverage.read_variable_sea_surface_wave_peak_period_at_time(
                    time_index)

                time_index += 1

        else:
            raise CoverageError("DefaultWriter","The given coverage is not an instance of 'TimeCoverage' or 'TimeLevelCoverage'")

    def write_variable_sea_surface_wave_from_direction(self):

        if (isinstance(self.coverage, TimeCoverage) or isinstance(self.coverage, TimeLevelCoverage)):

            var = self.ncfile.createVariable(VariableDefinition.VARIABLE_NAME['sea_surface_wave_from_direction'], float32,
                                             (VariableDefinition.VARIABLE_NAME['time'], VariableDefinition.VARIABLE_NAME['latitude'], VariableDefinition.VARIABLE_NAME['longitude'],), fill_value=9.96921e+36)
            var.long_name = VariableDefinition.LONG_NAME['sea_surface_wave_from_direction']
            var.standard_name = VariableDefinition.STANDARD_NAME['sea_surface_wave_from_direction']
            var.units = VariableDefinition.CANONICAL_UNITS['sea_surface_wave_from_direction']

            if self.coverage.rank == 0:
                logging.info('[DefaultWriter] Writing variable \'' + str(VariableDefinition.LONG_NAME['sea_surface_wave_from_direction']) + '\'')

            time_index = 0
            for time in self.coverage.read_axis_t():
                logging.debug(
                    '[DefaultWriter] Writing variable \'' + str(VariableDefinition.LONG_NAME['sea_surface_wave_from_direction']) + '\' at time \'' + str(time) + '\'')

                var[
                self.coverage.map_mpi[self.coverage.rank]["dst_global_t"].start+time_index:self.coverage.map_mpi[self.coverage.rank]["dst_global_t"].start+time_index+1,
                self.coverage.map_mpi[self.coverage.rank]["dst_global_y"],
                self.coverage.map_mpi[self.coverage.rank]["dst_global_x"]
                ] = self.coverage.read_variable_sea_surface_wave_from_direction_at_time(
                    time_index)

                time_index += 1

        else:
            raise CoverageError("DefaultWriter","The given coverage is not an instance of 'TimeCoverage' or 'TimeLevelCoverage'")

    def write_variable_sea_surface_wave_to_direction(self):

        if (isinstance(self.coverage, TimeCoverage) or isinstance(self.coverage, TimeLevelCoverage)):

            var = self.ncfile.createVariable(VariableDefinition.VARIABLE_NAME['sea_surface_wave_to_direction'], float32,
                                             (VariableDefinition.VARIABLE_NAME['time'], VariableDefinition.VARIABLE_NAME['latitude'], VariableDefinition.VARIABLE_NAME['longitude'],), fill_value=9.96921e+36)
            var.long_name = VariableDefinition.LONG_NAME['sea_surface_wave_to_direction']
            var.standard_name = VariableDefinition.STANDARD_NAME['sea_surface_wave_to_direction']
            var.units =VariableDefinition.CANONICAL_UNITS['sea_surface_wave_to_direction']

            if self.coverage.rank == 0:
                logging.info('[DefaultWriter] Writing variable \'' + str(VariableDefinition.LONG_NAME['sea_surface_wave_to_direction']) + '\'')

            time_index = 0
            for time in self.coverage.read_axis_t():
                logging.debug(
                    '[DefaultWriter] Writing variable \'' + str(VariableDefinition.LONG_NAME['sea_surface_wave_to_direction']) + '\' at time \'' + str(time) + '\'')

                var[
                self.coverage.map_mpi[self.coverage.rank]["dst_global_t"].start+time_index:self.coverage.map_mpi[self.coverage.rank]["dst_global_t"].start+time_index+1,
                self.coverage.map_mpi[self.coverage.rank]["dst_global_y"],
                self.coverage.map_mpi[self.coverage.rank]["dst_global_x"]
                ] = self.coverage.read_variable_sea_surface_wave_to_direction_at_time(
                    time_index)

                time_index += 1

        else:
            raise CoverageError("DefaultWriter","The given coverage is not an instance of 'TimeCoverage' or 'TimeLevelCoverage'")

    def write_variable_sea_surface_wave_stokes_drift_velocity(self):

        if (isinstance(self.coverage, TimeCoverage) or isinstance(self.coverage, TimeLevelCoverage)):

            ucur = self.ncfile.createVariable(VariableDefinition.VARIABLE_NAME['eastward_sea_surface_wave_stokes_drift_velocity'], float32, (VariableDefinition.VARIABLE_NAME['time'],VariableDefinition.VARIABLE_NAME['latitude'], VariableDefinition.VARIABLE_NAME['longitude'],),fill_value=9.96921e+36)
            ucur.long_name = VariableDefinition.LONG_NAME['eastward_sea_surface_wave_stokes_drift_velocity']
            ucur.standard_name = VariableDefinition.STANDARD_NAME['eastward_sea_surface_wave_stokes_drift_velocity']
            ucur.units = VariableDefinition.CANONICAL_UNITS['eastward_sea_surface_wave_stokes_drift_velocity']

            vcur = self.ncfile.createVariable(VariableDefinition.VARIABLE_NAME['northward_sea_surface_wave_stokes_drift_velocity'], float32, (VariableDefinition.VARIABLE_NAME['time'], VariableDefinition.VARIABLE_NAME['latitude'], VariableDefinition.VARIABLE_NAME['longitude'],),fill_value=9.96921e+36)
            vcur.long_name = VariableDefinition.LONG_NAME['northward_sea_surface_wave_stokes_drift_velocity']
            vcur.standard_name = VariableDefinition.STANDARD_NAME['northward_sea_surface_wave_stokes_drift_velocity']
            vcur.units = VariableDefinition.CANONICAL_UNITS['northward_sea_surface_wave_stokes_drift_velocity']

            if self.coverage.rank == 0:
                logging.info('[DefaultWriter] Writing variable \'Surface Stokes Drift Velocity\'\'')

            time_index=0
            for time in self.coverage.read_axis_t():

                logging.debug('[DefaultWriter] Writing variable \'Surface Stokes Drift Velocity\' at time \''+str(time)+'\'')

                cur = self.coverage.read_variable_sea_surface_wave_stokes_drift_velocity_at_time(time_index)

                ucur[
                self.coverage.map_mpi[self.coverage.rank]["dst_global_t"].start+time_index:self.coverage.map_mpi[self.coverage.rank]["dst_global_t"].start+time_index+1,
                self.coverage.map_mpi[self.coverage.rank]["dst_global_y"],
                self.coverage.map_mpi[self.coverage.rank]["dst_global_x"]
                ] = cur[0]

                vcur[
                self.coverage.map_mpi[self.coverage.rank]["dst_global_t"].start+time_index:self.coverage.map_mpi[self.coverage.rank]["dst_global_t"].start+time_index+1,
                self.coverage.map_mpi[self.coverage.rank]["dst_global_y"],
                self.coverage.map_mpi[self.coverage.rank]["dst_global_x"]
                ] = cur[1]

                time_index += 1

        else:
            raise CoverageError("DefaultWriter","The given coverage is not an instance of 'TimeCoverage' or 'TimeLevelCoverage'")

    def write_variable_radiation_pressure_bernouilli_head(self):

        if (isinstance(self.coverage, TimeCoverage) or isinstance(self.coverage, TimeLevelCoverage)):

            var = self.ncfile.createVariable(VariableDefinition.VARIABLE_NAME['radiation_pressure_bernouilli_head'], float32,
                                             (VariableDefinition.VARIABLE_NAME['time'], VariableDefinition.VARIABLE_NAME['latitude'], VariableDefinition.VARIABLE_NAME['longitude'],), fill_value=9.96921e+36)
            var.long_name = VariableDefinition.LONG_NAME['radiation_pressure_bernouilli_head']
            var.standard_name = VariableDefinition.STANDARD_NAME['radiation_pressure_bernouilli_head']
            var.units = VariableDefinition.CANONICAL_UNITS['radiation_pressure_bernouilli_head']

            if self.coverage.rank == 0:
                logging.info('[DefaultWriter] Writing variable \'' + str(VariableDefinition.LONG_NAME['radiation_pressure_bernouilli_head']) + '\'')

            time_index = 0
            for time in self.coverage.read_axis_t():
                logging.debug(
                    '[DefaultWriter] Writing variable \'' + str(VariableDefinition.LONG_NAME['radiation_pressure_bernouilli_head']) + '\' at time \'' + str(time) + '\'')

                var[
                self.coverage.map_mpi[self.coverage.rank]["dst_global_t"].start+time_index:self.coverage.map_mpi[self.coverage.rank]["dst_global_t"].start+time_index+1,
                self.coverage.map_mpi[self.coverage.rank]["dst_global_y"],
                self.coverage.map_mpi[self.coverage.rank]["dst_global_x"]
                ] = self.coverage.read_variable_radiation_pressure_bernouilli_head_at_time(time_index)

                time_index += 1

        else:
            raise CoverageError("DefaultWriter","The given coverage is not an instance of 'TimeCoverage' or 'TimeLevelCoverage'")

    def write_variable_sea_surface_wave_energy_flux_to_ocean(self):

        if (isinstance(self.coverage, TimeCoverage) or isinstance(self.coverage, TimeLevelCoverage)):

            var = self.ncfile.createVariable(VariableDefinition.VARIABLE_NAME['sea_surface_wave_energy_flux_to_ocean'],
                                             float32,
                                             (VariableDefinition.VARIABLE_NAME['time'],
                                              VariableDefinition.VARIABLE_NAME['latitude'],
                                              VariableDefinition.VARIABLE_NAME['longitude'],), fill_value=9.96921e+36)
            var.long_name = VariableDefinition.LONG_NAME['sea_surface_wave_energy_flux_to_ocean']
            var.standard_name = VariableDefinition.STANDARD_NAME['sea_surface_wave_energy_flux_to_ocean']
            var.units = VariableDefinition.CANONICAL_UNITS['sea_surface_wave_energy_flux_to_ocean']

            if self.coverage.rank == 0:
                logging.info('[DefaultWriter] Writing variable \'' + str(
                VariableDefinition.LONG_NAME['sea_surface_wave_energy_flux_to_ocean']) + '\'')

            time_index = 0
            for time in self.coverage.read_axis_t():
                logging.debug(
                    '[DefaultWriter] Writing variable \'' + str(
                        VariableDefinition.LONG_NAME['sea_surface_wave_energy_flux_to_ocean']) + '\' at time \'' + str(
                        time) + '\'')

                var[
                self.coverage.map_mpi[self.coverage.rank]["dst_global_t"].start+time_index:self.coverage.map_mpi[self.coverage.rank]["dst_global_t"].start+time_index+1,
                self.coverage.map_mpi[self.coverage.rank]["dst_global_y"],
                self.coverage.map_mpi[self.coverage.rank]["dst_global_x"]
                ] = self.coverage.read_variable_sea_surface_wave_energy_flux_to_ocean_at_time(time_index)

                time_index += 1

        else:
            raise CoverageError("DefaultWriter","The given coverage is not an instance of 'TimeCoverage' or 'TimeLevelCoverage'")

    #################
    # WAVES
    # Ground level
    #################
    def write_variable_sea_surface_wave_energy_dissipation_at_ground_level(self):

        if (isinstance(self.coverage, TimeCoverage) or isinstance(self.coverage, TimeLevelCoverage)):

            var = self.ncfile.createVariable(
                VariableDefinition.VARIABLE_NAME['sea_surface_wave_energy_dissipation_at_ground_level'], float32,
                (VariableDefinition.VARIABLE_NAME['time'], VariableDefinition.VARIABLE_NAME['latitude'],
                 VariableDefinition.VARIABLE_NAME['longitude'],), fill_value=9.96921e+36)
            var.long_name = VariableDefinition.LONG_NAME['sea_surface_wave_energy_dissipation_at_ground_level']
            var.standard_name = VariableDefinition.STANDARD_NAME['sea_surface_wave_energy_dissipation_at_ground_level']
            var.units = VariableDefinition.CANONICAL_UNITS['sea_surface_wave_energy_dissipation_at_ground_level']

            if self.coverage.rank == 0:
                logging.info('[DefaultWriter] Writing variable \'' + str(
                VariableDefinition.LONG_NAME['sea_surface_wave_energy_dissipation_at_ground_level']) + '\'')

            time_index = 0
            for time in self.coverage.read_axis_t():
                logging.debug(
                    '[DefaultWriter] Writing variable \'' + str(VariableDefinition.LONG_NAME[
                                                                    'sea_surface_wave_energy_dissipation_at_ground_level']) + '\' at time \'' + str(
                        time) + '\'')

                var[
                self.coverage.map_mpi[self.coverage.rank]["dst_global_t"].start+time_index:self.coverage.map_mpi[self.coverage.rank]["dst_global_t"].start+time_index+1,
                self.coverage.map_mpi[self.coverage.rank]["dst_global_y"],
                self.coverage.map_mpi[self.coverage.rank]["dst_global_x"]
                ] = self.coverage.read_variable_sea_surface_wave_energy_dissipation_at_ground_level_at_time(time_index)

                time_index += 1

        else:
            raise CoverageError("DefaultWriter","The given coverage is not an instance of 'TimeCoverage' or 'TimeLevelCoverage'")

    #################
    # WAVES
    # Momentum flux
    #################
    def write_variable_atmosphere_momentum_flux_to_waves(self):

        if (isinstance(self.coverage, TimeCoverage) or isinstance(self.coverage, TimeLevelCoverage)):

            ucur = self.ncfile.createVariable(VariableDefinition.VARIABLE_NAME['eastward_atmosphere_momentum_flux_to_waves'], float32, (VariableDefinition.VARIABLE_NAME['time'],VariableDefinition.VARIABLE_NAME['latitude'], VariableDefinition.VARIABLE_NAME['longitude'],),fill_value=9.96921e+36)
            ucur.long_name = VariableDefinition.LONG_NAME['eastward_atmosphere_momentum_flux_to_waves']
            ucur.standard_name = VariableDefinition.STANDARD_NAME['eastward_atmosphere_momentum_flux_to_waves']
            ucur.units = VariableDefinition.CANONICAL_UNITS['eastward_atmosphere_momentum_flux_to_waves']

            vcur = self.ncfile.createVariable(VariableDefinition.VARIABLE_NAME['northward_atmosphere_momentum_flux_to_waves'], float32, (VariableDefinition.VARIABLE_NAME['time'], VariableDefinition.VARIABLE_NAME['latitude'], VariableDefinition.VARIABLE_NAME['longitude'],),fill_value=9.96921e+36)
            vcur.long_name = VariableDefinition.LONG_NAME['northward_atmosphere_momentum_flux_to_waves']
            vcur.standard_name = VariableDefinition.STANDARD_NAME['northward_atmosphere_momentum_flux_to_waves']
            vcur.units =VariableDefinition.CANONICAL_UNITS['northward_atmosphere_momentum_flux_to_waves']

            if self.coverage.rank == 0:
                logging.info('[DefaultWriter] Writing variable \'Atmosphere Momentum Flux to Waves\'\'')

            time_index=0
            for time in self.coverage.read_axis_t():

                logging.debug('[DefaultWriter] Writing variable \'Atmosphere Momentum Flux to Waves\' at time \''+str(time)+'\'')

                cur = self.coverage.read_variable_atmosphere_momentum_flux_to_waves_at_time(time_index)

                ucur[
                self.coverage.map_mpi[self.coverage.rank]["dst_global_t"].start+time_index:self.coverage.map_mpi[self.coverage.rank]["dst_global_t"].start+time_index+1,
                self.coverage.map_mpi[self.coverage.rank]["dst_global_y"],
                self.coverage.map_mpi[self.coverage.rank]["dst_global_x"]
                ] = cur[0]

                vcur[
                self.coverage.map_mpi[self.coverage.rank]["dst_global_t"].start+time_index:self.coverage.map_mpi[self.coverage.rank]["dst_global_t"].start+time_index+1,
                self.coverage.map_mpi[self.coverage.rank]["dst_global_y"],
                self.coverage.map_mpi[self.coverage.rank]["dst_global_x"]
                ] = cur[1]

                time_index += 1

        else:
            raise CoverageError("DefaultWriter","The given coverage is not an instance of 'TimeCoverage' or 'TimeLevelCoverage'")

    def write_variable_waves_momentum_flux_to_ocean(self):

        if (isinstance(self.coverage, TimeCoverage) or isinstance(self.coverage, TimeLevelCoverage)):

            ucur = self.ncfile.createVariable(VariableDefinition.VARIABLE_NAME['eastward_waves_momentum_flux_to_ocean'], float32, (VariableDefinition.VARIABLE_NAME['time'],VariableDefinition.VARIABLE_NAME['latitude'], VariableDefinition.VARIABLE_NAME['longitude'],),fill_value=9.96921e+36)
            ucur.long_name = VariableDefinition.LONG_NAME['eastward_waves_momentum_flux_to_ocean']
            ucur.standard_name = VariableDefinition.STANDARD_NAME['eastward_waves_momentum_flux_to_ocean']
            ucur.units = VariableDefinition.CANONICAL_UNITS['eastward_waves_momentum_flux_to_ocean']

            vcur = self.ncfile.createVariable(VariableDefinition.VARIABLE_NAME['northward_waves_momentum_flux_to_ocean'], float32, (VariableDefinition.VARIABLE_NAME['time'], VariableDefinition.VARIABLE_NAME['latitude'], VariableDefinition.VARIABLE_NAME['longitude'],),fill_value=9.96921e+36)
            vcur.long_name = VariableDefinition.LONG_NAME['northward_waves_momentum_flux_to_ocean']
            vcur.standard_name = VariableDefinition.STANDARD_NAME['northward_waves_momentum_flux_to_ocean']
            vcur.units = VariableDefinition.CANONICAL_UNITS['northward_waves_momentum_flux_to_ocean']

            if self.coverage.rank == 0:
                logging.info('[DefaultWriter] Writing variable \'Waves Momentum Flux To Ocean\'\'')

            time_index=0
            for time in self.coverage.read_axis_t():

                logging.debug('[DefaultWriter] Writing variable \'Waves Momentum Flux To Ocean\' at time \''+str(time)+'\'')

                cur = self.coverage.read_variable_waves_momentum_flux_to_ocean_at_time(time_index)

                ucur[
                self.coverage.map_mpi[self.coverage.rank]["dst_global_t"].start+time_index:self.coverage.map_mpi[self.coverage.rank]["dst_global_t"].start+time_index+1,
                self.coverage.map_mpi[self.coverage.rank]["dst_global_y"],
                self.coverage.map_mpi[self.coverage.rank]["dst_global_x"]
                ] = cur[0]

                vcur[
                self.coverage.map_mpi[self.coverage.rank]["dst_global_t"].start+time_index:self.coverage.map_mpi[self.coverage.rank]["dst_global_t"].start+time_index+1,
                self.coverage.map_mpi[self.coverage.rank]["dst_global_y"],
                self.coverage.map_mpi[self.coverage.rank]["dst_global_x"]
                ] = cur[1]

                time_index += 1

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
    def write_variable_wind_stress(self):

        if (isinstance(self.coverage, TimeCoverage) or isinstance(self.coverage, TimeLevelCoverage)):

            ucur = self.ncfile.createVariable(VariableDefinition.VARIABLE_NAME['eastward_wind_stress'], float32, (
            VariableDefinition.VARIABLE_NAME['time'], VariableDefinition.VARIABLE_NAME['latitude'],
            VariableDefinition.VARIABLE_NAME['longitude'],),
                                              fill_value=9.96921e+36)
            ucur.long_name = VariableDefinition.LONG_NAME['eastward_wind_stress']
            ucur.standard_name = VariableDefinition.STANDARD_NAME['eastward_wind_stress']
            ucur.units = VariableDefinition.CANONICAL_UNITS['eastward_wind_stress']

            vcur = self.ncfile.createVariable(VariableDefinition.VARIABLE_NAME['northward_wind_stress'], float32, (
            VariableDefinition.VARIABLE_NAME['time'], VariableDefinition.VARIABLE_NAME['latitude'],
            VariableDefinition.VARIABLE_NAME['longitude'],),
                                              fill_value=9.96921e+36)
            vcur.long_name = VariableDefinition.LONG_NAME['northward_wind_stress']
            vcur.standard_name = VariableDefinition.STANDARD_NAME['northward_wind_stress']
            vcur.units = VariableDefinition.CANONICAL_UNITS['northward_wind_stress']

            if self.coverage.rank == 0:
                logging.info('[DefaultWriter] Writing variable \'Wind Stress\'\'')

            time_index = 0
            for time in self.coverage.read_axis_t():
                logging.debug('[DefaultWriter] Writing variable \'Wind Stress\' at time \'' + str(time) + '\'')

                cur = self.coverage.read_variable_wind_stress_at_time(time_index)

                ucur[
                self.coverage.map_mpi[self.coverage.rank]["dst_global_t"].start+time_index:self.coverage.map_mpi[self.coverage.rank]["dst_global_t"].start+time_index+1,
                self.coverage.map_mpi[self.coverage.rank]["dst_global_y"],
                self.coverage.map_mpi[self.coverage.rank]["dst_global_x"]
                ] = cur[0]

                vcur[
                self.coverage.map_mpi[self.coverage.rank]["dst_global_t"].start+time_index:self.coverage.map_mpi[self.coverage.rank]["dst_global_t"].start+time_index+1,
                self.coverage.map_mpi[self.coverage.rank]["dst_global_y"],
                self.coverage.map_mpi[self.coverage.rank]["dst_global_x"]
                ] = cur[1]

                time_index += 1

        else:
            raise CoverageError("DefaultWriter","The given coverage is not an instance of 'TimeCoverage' or 'TimeLevelCoverage'")

    #################
    # METEO
    # At 10 m
    #################
    def write_variable_wind_10m(self):

        if (isinstance(self.coverage, TimeCoverage) or isinstance(self.coverage, TimeLevelCoverage)):

            ucur = self.ncfile.createVariable(VariableDefinition.VARIABLE_NAME['eastward_wind_10m'], float32, (VariableDefinition.VARIABLE_NAME['time'],VariableDefinition.VARIABLE_NAME['latitude'], VariableDefinition.VARIABLE_NAME['longitude'],),fill_value=9.96921e+36)
            ucur.long_name = VariableDefinition.LONG_NAME['eastward_wind_10m']
            ucur.standard_name = VariableDefinition.STANDARD_NAME['eastward_wind_10m']
            ucur.units = VariableDefinition.CANONICAL_UNITS['eastward_wind_10m']

            vcur = self.ncfile.createVariable(VariableDefinition.VARIABLE_NAME['northward_wind_10m'], float32, (VariableDefinition.VARIABLE_NAME['time'], VariableDefinition.VARIABLE_NAME['latitude'], VariableDefinition.VARIABLE_NAME['longitude'],),fill_value=9.96921e+36)
            vcur.long_name = VariableDefinition.LONG_NAME['northward_wind_10m']
            vcur.standard_name = VariableDefinition.STANDARD_NAME['northward_wind_10m']
            vcur.units = VariableDefinition.CANONICAL_UNITS['northward_wind_10m']

            if self.coverage.rank == 0:
                logging.info('[DefaultWriter] Writing variable \'Wind 10m\'\'')

            time_index=0
            for time in self.coverage.read_axis_t():

                logging.debug('[DefaultWriter] Writing variable \'Wind 10m\' at time \''+str(time)+'\'')

                cur = self.coverage.read_variable_wind_10m_at_time(time_index)

                ucur[
                self.coverage.map_mpi[self.coverage.rank]["dst_global_t"].start+time_index:self.coverage.map_mpi[self.coverage.rank]["dst_global_t"].start+time_index+1,
                self.coverage.map_mpi[self.coverage.rank]["dst_global_y"],
                self.coverage.map_mpi[self.coverage.rank]["dst_global_x"]
                ] = cur[0]

                vcur[
                self.coverage.map_mpi[self.coverage.rank]["dst_global_t"].start+time_index:self.coverage.map_mpi[self.coverage.rank]["dst_global_t"].start+time_index+1,
                self.coverage.map_mpi[self.coverage.rank]["dst_global_y"],
                self.coverage.map_mpi[self.coverage.rank]["dst_global_x"]
                ] = cur[1]

                time_index += 1

        else:
            raise CoverageError("DefaultWriter","The given coverage is not an instance of 'TimeCoverage' or 'TimeLevelCoverage'")

    def write_variable_wind_speed_10m(self):

        if (isinstance(self.coverage, TimeCoverage) or isinstance(self.coverage, TimeLevelCoverage)):

            var = self.ncfile.createVariable(VariableDefinition.VARIABLE_NAME['wind_speed_10m'], float32,
                                             (VariableDefinition.VARIABLE_NAME['time'], VariableDefinition.VARIABLE_NAME['latitude'], VariableDefinition.VARIABLE_NAME['longitude'],), fill_value=9.96921e+36)
            var.long_name = VariableDefinition.LONG_NAME['wind_speed_10m']
            var.standard_name = VariableDefinition.STANDARD_NAME['wind_speed_10m']
            var.units = VariableDefinition.CANONICAL_UNITS['wind_speed_10m']

            if self.coverage.rank == 0:
                logging.info('[DefaultWriter] Writing variable \'' + str(VariableDefinition.LONG_NAME['wind_speed_10m']) + '\'')

            time_index = 0
            for time in self.coverage.read_axis_t():
                logging.debug(
                    '[DefaultWriter] Writing variable \'' + str(VariableDefinition.LONG_NAME['wind_speed_10m']) + '\' at time \'' + str(time) + '\'')

                var[
                self.coverage.map_mpi[self.coverage.rank]["dst_global_t"].start+time_index:self.coverage.map_mpi[self.coverage.rank]["dst_global_t"].start+time_index+1,
                self.coverage.map_mpi[self.coverage.rank]["dst_global_y"],
                self.coverage.map_mpi[self.coverage.rank]["dst_global_x"]
                ] = self.coverage.read_variable_wind_speed_10m_at_time(time_index)

                time_index += 1

        else:
            raise CoverageError("DefaultWriter","The given coverage is not an instance of 'TimeCoverage' or 'TimeLevelCoverage'")

    def write_variable_wind_to_direction_10m(self):

        if (isinstance(self.coverage, TimeCoverage) or isinstance(self.coverage, TimeLevelCoverage)):

            var = self.ncfile.createVariable(VariableDefinition.VARIABLE_NAME['wind_to_direction_10m'], float32,
                                             (VariableDefinition.VARIABLE_NAME['time'], VariableDefinition.VARIABLE_NAME['latitude'], VariableDefinition.VARIABLE_NAME['longitude'],), fill_value=9.96921e+36)
            var.long_name = VariableDefinition.LONG_NAME['wind_to_direction_10m']
            var.standard_name = VariableDefinition.STANDARD_NAME['wind_to_direction_10m']
            var.units = VariableDefinition.CANONICAL_UNITS['wind_to_direction_10m']

            if self.coverage.rank == 0:
                logging.info('[DefaultWriter] Writing variable \'' + str(VariableDefinition.LONG_NAME['wind_to_direction_10m']) + '\'')

            time_index = 0
            for time in self.coverage.read_axis_t():
                logging.debug(
                    '[DefaultWriter] Writing variable \'' + str(VariableDefinition.LONG_NAME['wind_to_direction_10m']) + '\' at time \'' + str(time) + '\'')

                var[
                self.coverage.map_mpi[self.coverage.rank]["dst_global_t"].start+time_index:self.coverage.map_mpi[self.coverage.rank]["dst_global_t"].start+time_index+1,
                self.coverage.map_mpi[self.coverage.rank]["dst_global_y"],
                self.coverage.map_mpi[self.coverage.rank]["dst_global_x"]
                ] = self.coverage.read_variable_wind_to_direction_10m_at_time(time_index)

                time_index += 1

        else:
            raise CoverageError("DefaultWriter","The given coverage is not an instance of 'TimeCoverage' or 'TimeLevelCoverage'")

    def write_variable_wind_from_direction_10m(self):

        if (isinstance(self.coverage, TimeCoverage) or isinstance(self.coverage, TimeLevelCoverage)):

            var = self.ncfile.createVariable(VariableDefinition.VARIABLE_NAME['wind_from_direction_10m'], float32,
                                             (VariableDefinition.VARIABLE_NAME['time'], VariableDefinition.VARIABLE_NAME['latitude'], VariableDefinition.VARIABLE_NAME['longitude'],), fill_value=9.96921e+36)
            var.long_name = VariableDefinition.LONG_NAME['wind_from_direction_10m']
            var.standard_name = VariableDefinition.STANDARD_NAME['wind_from_direction_10m']
            var.units = VariableDefinition.CANONICAL_UNITS['wind_from_direction_10m']

            if self.coverage.rank == 0:
                logging.info('[DefaultWriter] Writing variable \'' + str(VariableDefinition.LONG_NAME['wind_from_direction_10m']) + '\'')

            time_index = 0
            for time in self.coverage.read_axis_t():
                logging.debug(
                    '[DefaultWriter] Writing variable \'' + str(VariableDefinition.LONG_NAME['wind_from_direction_10m']) + '\' at time \'' + str(time) + '\'')

                var[
                self.coverage.map_mpi[self.coverage.rank]["dst_global_t"].start+time_index:self.coverage.map_mpi[self.coverage.rank]["dst_global_t"].start+time_index+1,
                self.coverage.map_mpi[self.coverage.rank]["dst_global_y"],
                self.coverage.map_mpi[self.coverage.rank]["dst_global_x"]
                ] = self.coverage.read_variable_wind_from_direction_10m_at_time(time_index)

                time_index += 1

        else:
            raise CoverageError("DefaultWriter","The given coverage is not an instance of 'TimeCoverage' or 'TimeLevelCoverage'")





