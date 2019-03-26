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
from numpy import float32
from numpy import float64
import numpy as np
import logging
from coverage.TimeCoverage import TimeCoverage
from coverage.LevelCoverage import LevelCoverage
from coverage.TimeLevelCoverage import TimeLevelCoverage

class DefaultWriter (CoverageWriter):

    def __init__(self,cov,myFile):
        CoverageWriter.__init__(self,cov,myFile);
        self.ncfile = Dataset(self.filename, 'w', format='NETCDF4')
        self.ncfile.description = 'Generated with pyGeoSpatialETL'

        # dimensions
        self.ncfile.createDimension(VariableDefinition.VARIABLE_NAME['latitude'], self.coverage.get_y_size())
        self.ncfile.createDimension(VariableDefinition.VARIABLE_NAME['longitude'], self.coverage.get_x_size())

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

            # data
            latitudes[:] = self.coverage.read_axis_y();
            longitudes[:] = self.coverage.read_axis_x();

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
            latitudes[:,:] = self.coverage.read_axis_y();
            longitudes[:,:] = self.coverage.read_axis_x();


        if(isinstance(self.coverage, TimeCoverage) or isinstance(self.coverage, TimeLevelCoverage)):

            self.ncfile.createDimension(VariableDefinition.VARIABLE_NAME['time'], None)
            times = self.ncfile.createVariable(VariableDefinition.VARIABLE_NAME['time'], float64, (VariableDefinition.VARIABLE_NAME['time'],))
            times.units= 'seconds since 1970-01-01 00:00:00'
            times.calendar= 'gregorian'
            times.standard_name= 'time'
            times.axis='T'
            times.conventions = "UTC time"

            times[:] = date2num(self.coverage.read_axis_t(), units = times.units, calendar = times.calendar)

        if(isinstance(self.coverage, LevelCoverage) or isinstance(self.coverage, TimeLevelCoverage)):

            if self.coverage.is_sigma_coordinate():
                raise ValueError("This writer supports only Coverage with a regular vertical axis.")

            self.ncfile.createDimension(VariableDefinition.VARIABLE_NAME['depth'], self.coverage.get_z_size())
            levels = self.ncfile.createVariable(VariableDefinition.VARIABLE_NAME['depth'], float64, (VariableDefinition.VARIABLE_NAME['depth'],))
            levels.standard_name= 'depth'
            levels.long_name="Positive depth"
            levels.axis='Z'
            levels.units = VariableDefinition.CANONICAL_UNITS['depth'];

            levels[:] = self.coverage.read_axis_z()

        
    def close(self):
        self.ncfile.close()

    # HYDRO
    def write_variable_mesh_size(self):

        var = self.ncfile.createVariable(VariableDefinition.VARIABLE_NAME['mesh_size'], float32, (VariableDefinition.VARIABLE_NAME['latitude'], VariableDefinition.VARIABLE_NAME['longitude'],), fill_value=9.96921e+36)
        var.long_name = VariableDefinition.LONG_NAME['mesh_size']
        var.standard_name = VariableDefinition.STANDARD_NAME['mesh_size']
        var.units = VariableDefinition.CANONICAL_UNITS['mesh_size']

        logging.info('[DefaultWriter] Writing variable \'' + str(VariableDefinition.LONG_NAME['mesh_size']) + '\'')

        var[:] = self.coverage.read_variable_mesh_size()

    def write_variable_2D_sea_binary_mask(self):

        var = self.ncfile.createVariable(VariableDefinition.VARIABLE_NAME['2d_sea_binary_mask'], int, (VariableDefinition.VARIABLE_NAME['latitude'],VariableDefinition.VARIABLE_NAME['longitude'],), fill_value=-9999)
        var.long_name = VariableDefinition.LONG_NAME['2d_sea_binary_mask']
        var.standard_name = VariableDefinition.STANDARD_NAME['2d_sea_binary_mask']
        var.units = VariableDefinition.CANONICAL_UNITS['2d_sea_binary_mask']
        var.comment = "1 = sea, 0 = land"

        logging.info('[DefaultWriter] Writing variable \'' + str(VariableDefinition.LONG_NAME['2d_sea_binary_mask']) + '\'')
        var[:] = self.coverage.read_variable_2D_sea_binary_mask()

    def write_variable_wet_binary_mask(self):

        var = self.ncfile.createVariable(VariableDefinition.VARIABLE_NAME['wet_binary_mask'], float32, (VariableDefinition.VARIABLE_NAME['time'], VariableDefinition.VARIABLE_NAME['latitude'], VariableDefinition.VARIABLE_NAME['longitude'],),
                                         fill_value=9.96921e+36)
        var.long_name = VariableDefinition.LONG_NAME['wet_binary_mask']
        var.standard_name = VariableDefinition.STANDARD_NAME['wet_binary_mask']
        var.units = VariableDefinition.CANONICAL_UNITS['wet_binary_mask']
        var.comment = "1 = sea, 0 = land"

        logging.info('[DefaultWriter] Writing variable \'' + str(VariableDefinition.VARIABLE_NAME['wet_binary_mask']) + '\'')

        time_index = 0
        for time in self.coverage.read_axis_t():
            logging.debug(
                '[DefaultWriter] Writing variable \'' + str(VariableDefinition.VARIABLE_NAME['wet_binary_mask']) + '\' at time \'' + str(time) + '\'')
            var[time_index:time_index + 1, :] = self.coverage.read_variable_2D_wet_binary_mask_at_time(time)
            time_index += 1

    def write_variable_bathymetry(self):

        var = self.ncfile.createVariable(VariableDefinition.VARIABLE_NAME['bathymetry'], float32, (VariableDefinition.VARIABLE_NAME['latitude'], VariableDefinition.VARIABLE_NAME['longitude'],),fill_value=9.96921e+36)
        var.long_name = VariableDefinition.LONG_NAME['bathymetry']
        var.standard_name = VariableDefinition.STANDARD_NAME['bathymetry']
        var.units = VariableDefinition.CANONICAL_UNITS['bathymetry']

        logging.info('[DefaultWriter] Writing variable \'' + str(VariableDefinition.LONG_NAME['bathymetry']) + '\'')

        var[:] = self.coverage.read_variable_bathymetry()
            
    def write_variable_sea_surface_height_above_mean_sea_level(self):
            
        var = self.ncfile.createVariable(VariableDefinition.VARIABLE_NAME['sea_surface_height_above_mean_sea_level'], float32, (VariableDefinition.VARIABLE_NAME['time'], VariableDefinition.VARIABLE_NAME['latitude'], VariableDefinition.VARIABLE_NAME['longitude'],),fill_value=9.96921e+36)
        var.long_name = VariableDefinition.LONG_NAME['sea_surface_height_above_mean_sea_level']
        var.standard_name = VariableDefinition.STANDARD_NAME['sea_surface_height_above_mean_sea_level']
        var.units = VariableDefinition.CANONICAL_UNITS['sea_surface_height_above_mean_sea_level']

        logging.info('[DefaultWriter] Writing variable \''+str(VariableDefinition.LONG_NAME['sea_surface_height_above_mean_sea_level'])+'\'')
        
        time_index=0
        for time in self.coverage.read_axis_t():
            logging.debug('[DefaultWriter] Writing variable \''+str(VariableDefinition.LONG_NAME['sea_surface_height_above_mean_sea_level'])+'\' at time \''+str(time)+'\'')
            var[time_index:time_index+1,:] = self.coverage.read_variable_sea_surface_height_above_mean_sea_level_at_time(time)
            time_index += 1

    def write_variable_sea_surface_height_above_geoid(self):

        var = self.ncfile.createVariable(VariableDefinition.VARIABLE_NAME['sea_surface_height_above_geoid'], float32,
                                         (VariableDefinition.VARIABLE_NAME['time'], VariableDefinition.VARIABLE_NAME['latitude'], VariableDefinition.VARIABLE_NAME['longitude'],), fill_value=9.96921e+36)
        var.long_name = VariableDefinition.LONG_NAME['sea_surface_height_above_geoid']
        var.standard_name = VariableDefinition.STANDARD_NAME['sea_surface_height_above_geoid']
        var.units = VariableDefinition.CANONICAL_UNITS['sea_surface_height_above_geoid']

        logging.info('[DefaultWriter] Writing variable \'' + str(VariableDefinition.LONG_NAME['sea_surface_height_above_geoid']) + '\'')

        time_index = 0
        for time in self.coverage.read_axis_t():
            logging.debug(
                '[DefaultWriter] Writing variable \'' + str(VariableDefinition.LONG_NAME['sea_surface_height_above_geoid']) + '\' at time \'' + str(time) + '\'')
            var[time_index:time_index + 1,
            :] = self.coverage.read_variable_sea_surface_height_above_geoid_at_time(time)
            time_index += 1

    def write_variable_baroclinic_sea_water_velocity(self):

        ucur = self.ncfile.createVariable(VariableDefinition.VARIABLE_NAME['baroclinic_eastward_sea_water_velocity'], float32,
                                          (VariableDefinition.VARIABLE_NAME['time'],VariableDefinition.VARIABLE_NAME['depth'], VariableDefinition.VARIABLE_NAME['latitude'], VariableDefinition.VARIABLE_NAME['longitude'],), fill_value=9.96921e+36)
        ucur.long_name = VariableDefinition.LONG_NAME['baroclinic_eastward_sea_water_velocity']
        ucur.standard_name = VariableDefinition.STANDARD_NAME['baroclinic_eastward_sea_water_velocity']
        ucur.units = VariableDefinition.CANONICAL_UNITS['baroclinic_eastward_sea_water_velocity']
        ucur.comment = "cur=sqrt(U**2+V**2)";

        vcur = self.ncfile.createVariable(VariableDefinition.VARIABLE_NAME['baroclinic_northward_sea_water_velocity'], float32,
                                          (VariableDefinition.VARIABLE_NAME['time'],VariableDefinition.VARIABLE_NAME['depth'], VariableDefinition.VARIABLE_NAME['latitude'], VariableDefinition.VARIABLE_NAME['longitude'],), fill_value=9.96921e+36)
        vcur.long_name = VariableDefinition.LONG_NAME['baroclinic_northward_sea_water_velocity']
        vcur.standard_name = VariableDefinition.STANDARD_NAME['baroclinic_northward_sea_water_velocity']
        vcur.units = VariableDefinition.CANONICAL_UNITS['baroclinic_northward_sea_water_velocity']
        vcur.comment = "cur=sqrt(U**2+V**2)";

        logging.info('[DefaultWriter] Writing variable \'Baroclinic Sea Water Velocity\'\'')

        time_index = 0
        for time in self.coverage.read_axis_t():

            logging.debug(
                '[DefaultWriter] Writing variable \'Baroclinic Sea Water Velocity\' at time \'' + str(time) + '\'')

            level_index = 0
            for level in self.coverage.read_axis_z():
                cur = self.coverage.read_variable_baroclinic_sea_water_velocity_at_time_and_depth(time, level)

                ucur[time_index:time_index + 1, level_index:level_index + 1, :, :] = cur[0]
                vcur[time_index:time_index + 1, level_index:level_index + 1, :, :] = cur[1]
                level_index += 1

            time_index += 1

    def write_variable_barotropic_sea_water_velocity(self):
        ucur = self.ncfile.createVariable(VariableDefinition.VARIABLE_NAME['barotropic_eastward_sea_water_velocity'], float32,
                                          (VariableDefinition.VARIABLE_NAME['time'], VariableDefinition.VARIABLE_NAME['latitude'], VariableDefinition.VARIABLE_NAME['longitude'],), fill_value=9.96921e+36)
        ucur.long_name = VariableDefinition.LONG_NAME['barotropic_eastward_sea_water_velocity']
        ucur.standard_name = VariableDefinition.STANDARD_NAME['barotropic_eastward_sea_water_velocity']
        ucur.units = VariableDefinition.CANONICAL_UNITS['barotropic_eastward_sea_water_velocity']
        ucur.comment = "cur=sqrt(U**2+V**2)";

        vcur = self.ncfile.createVariable(VariableDefinition.VARIABLE_NAME['barotropic_northward_sea_water_velocity'], float32,
                                          (VariableDefinition.VARIABLE_NAME['time'], VariableDefinition.VARIABLE_NAME['latitude'], VariableDefinition.VARIABLE_NAME['longitude'],), fill_value=9.96921e+36)
        vcur.long_name = VariableDefinition.LONG_NAME['barotropic_northward_sea_water_velocity']
        vcur.standard_name = VariableDefinition.STANDARD_NAME['barotropic_northward_sea_water_velocity']
        vcur.units = VariableDefinition.CANONICAL_UNITS['barotropic_northward_sea_water_velocity']
        vcur.comment = "cur=sqrt(U**2+V**2)";

        logging.info('[DefaultWriter] Writing variable \'Barotropic Sea Water Velocity\'\'')

        time_index = 0
        for time in self.coverage.read_axis_t():
            logging.debug(
                '[DefaultWriter] Writing variable \'Barotropic Sea Water Velocity\' at time \'' + str(time) + '\'')

            cur = self.coverage.read_variable_barotropic_sea_water_velocity_at_time(time)

            ucur[time_index:time_index + 1, :, :] = cur[0]
            vcur[time_index:time_index + 1, :, :] = cur[1]
            time_index += 1

    def write_variable_sea_surface_temperature(self):

        var = self.ncfile.createVariable(VariableDefinition.VARIABLE_NAME['sea_surface_temperature'], float32, (VariableDefinition.VARIABLE_NAME['time'], VariableDefinition.VARIABLE_NAME['latitude'], VariableDefinition.VARIABLE_NAME['longitude'],),
                                         fill_value=9.96921e+36)
        var.long_name = VariableDefinition.LONG_NAME['sea_surface_temperature']
        var.standard_name = VariableDefinition.STANDARD_NAME['sea_surface_temperature']
        var.units = VariableDefinition.CANONICAL_UNITS['sea_surface_temperature']

        logging.info('[DefaultWriter] Writing variable \'' + str(VariableDefinition.LONG_NAME['sea_surface_temperature']) + '\'')

        time_index = 0
        for time in self.coverage.read_axis_t():
            logging.debug(
                '[DefaultWriter] Writing variable \'' + str(VariableDefinition.LONG_NAME['sea_surface_temperature']) + '\' at time \'' + str(time) + '\'')
            var[time_index:time_index + 1, :] = self.coverage.read_variable_sea_surface_temperature_at_time(time)
            time_index += 1

    def write_variable_sea_water_temperature(self):

        var = self.ncfile.createVariable(VariableDefinition.VARIABLE_NAME['sea_water_temperature'], float32,
                                         (VariableDefinition.VARIABLE_NAME['time'], VariableDefinition.VARIABLE_NAME['depth'], VariableDefinition.VARIABLE_NAME['latitude'], VariableDefinition.VARIABLE_NAME['longitude'],), fill_value=9.96921e+36)
        var.long_name = VariableDefinition.LONG_NAME['sea_water_temperature']
        var.standard_name = VariableDefinition.STANDARD_NAME['sea_water_temperature']
        var.units = VariableDefinition.CANONICAL_UNITS['sea_water_temperature']

        logging.info('[DefaultWriter] Writing variable \'' + str(VariableDefinition.LONG_NAME['sea_water_temperature']) + '\'')

        time_index = 0
        for time in self.coverage.read_axis_t():

            logging.debug(
                '[DefaultWriter] Writing variable \'' + str(VariableDefinition.LONG_NAME['sea_water_temperature']) + '\' at time \'' + str(time) + '\'')

            level_index = 0
            for level in self.coverage.read_axis_z():
                var[time_index:time_index + 1, level_index:level_index + 1, :,
                :] = self.coverage.read_variable_sea_water_temperature_at_time_and_depth(time, level)
                level_index += 1

            time_index += 1

    def write_variable_sea_water_salinity(self):

        var = self.ncfile.createVariable(VariableDefinition.VARIABLE_NAME['sea_water_salinity'], float32,
                                         (VariableDefinition.VARIABLE_NAME['time'], VariableDefinition.VARIABLE_NAME['depth'], VariableDefinition.VARIABLE_NAME['latitude'], VariableDefinition.VARIABLE_NAME['longitude'],), fill_value=9.96921e+36)
        var.long_name = VariableDefinition.LONG_NAME['sea_water_salinity']
        var.standard_name = VariableDefinition.STANDARD_NAME['sea_water_salinity']
        var.units = VariableDefinition.CANONICAL_UNITS['sea_water_salinity']

        logging.info('[DefaultWriter] Writing variable \'' + str(VariableDefinition.LONG_NAME['sea_water_salinity']) + '\'')

        time_index = 0
        for time in self.coverage.read_axis_t():

            logging.debug(
                '[DefaultWriter] Writing variable \'' + str(VariableDefinition.LONG_NAME['sea_water_salinity']) + '\' at time \'' + str(time) + '\'')

            level_index = 0
            for level in self.coverage.read_axis_z():
                var[time_index:time_index + 1, level_index:level_index + 1, :,
                :] = self.coverage.read_variable_sea_water_salinity_at_time_and_depth(time, level)
                level_index += 1

            time_index += 1

    # WAVES
    def write_variable_sea_surface_wave_significant_height(self):

        var = self.ncfile.createVariable(VariableDefinition.VARIABLE_NAME['sea_surface_wave_significant_height'], float32, (VariableDefinition.VARIABLE_NAME['time'], VariableDefinition.VARIABLE_NAME['latitude'], VariableDefinition.VARIABLE_NAME['longitude'],),fill_value=9.96921e+36)
        var.long_name = VariableDefinition.LONG_NAME['sea_surface_wave_significant_height']
        var.standard_name = VariableDefinition.STANDARD_NAME['sea_surface_wave_significant_height']
        var.units = VariableDefinition.CANONICAL_UNITS['sea_surface_wave_significant_height']

        logging.info('[DefaultWriter] Writing variable \'' + str(VariableDefinition.LONG_NAME['sea_surface_wave_significant_height']) + '\'')

        time_index=0
        for time in self.coverage.read_axis_t():
            logging.debug(
                '[DefaultWriter] Writing variable \'' + str(VariableDefinition.LONG_NAME['sea_surface_wave_significant_height']) + '\' at time \'' + str(time) + '\'')
            var[time_index:time_index+1,:] = self.coverage.read_variable_sea_surface_wave_significant_height_at_time(time)
            time_index += 1

    def write_variable_sea_surface_wave_breaking_height(self):

        var = self.ncfile.createVariable(VariableDefinition.VARIABLE_NAME['sea_surface_wave_breaking_height'], float32,
                                         (VariableDefinition.VARIABLE_NAME['time'], VariableDefinition.VARIABLE_NAME['latitude'], VariableDefinition.VARIABLE_NAME['longitude'],), fill_value=9.96921e+36)
        var.long_name = VariableDefinition.LONG_NAME['sea_surface_wave_breaking_height']
        var.standard_name = VariableDefinition.STANDARD_NAME['sea_surface_wave_breaking_height']
        var.units =VariableDefinition.CANONICAL_UNITS['sea_surface_wave_breaking_height']

        logging.info('[DefaultWriter] Writing variable \'' + str(VariableDefinition.LONG_NAME['sea_surface_wave_breaking_height']) + '\'')

        time_index = 0
        for time in self.coverage.read_axis_t():
            logging.debug(
                '[DefaultWriter] Writing variable \'' + str(VariableDefinition.LONG_NAME['sea_surface_wave_breaking_height']) + '\' at time \'' + str(time) + '\'')
            var[time_index:time_index + 1,
            :] = self.coverage.read_variable_sea_surface_wave_breaking_height_at_time(time)
            time_index += 1

    def write_variable_sea_surface_wave_mean_period(self):

        var = self.ncfile.createVariable(VariableDefinition.VARIABLE_NAME['sea_surface_wave_mean_period'], float32,
                                         (VariableDefinition.VARIABLE_NAME['time'], VariableDefinition.VARIABLE_NAME['latitude'], VariableDefinition.VARIABLE_NAME['longitude'],), fill_value=9.96921e+36)
        var.long_name = VariableDefinition.LONG_NAME['sea_surface_wave_mean_period']
        var.standard_name = VariableDefinition.STANDARD_NAME['sea_surface_wave_mean_period']
        var.units = VariableDefinition.CANONICAL_UNITS['sea_surface_wave_mean_period']

        logging.info('[DefaultWriter] Writing variable \'' + str(VariableDefinition.LONG_NAME['sea_surface_wave_mean_period']) + '\'')

        time_index = 0
        for time in self.coverage.read_axis_t():
            logging.debug(
                '[DefaultWriter] Writing variable \'' + str(VariableDefinition.LONG_NAME['sea_surface_wave_mean_period']) + '\' at time \'' + str(time) + '\'')
            var[time_index:time_index + 1, :] = self.coverage.read_variable_sea_surface_wave_mean_period_at_time(
                time)
            time_index += 1

    def write_variable_sea_surface_wave_peak_period(self):

        var = self.ncfile.createVariable(VariableDefinition.VARIABLE_NAME['sea_surface_wave_peak_period'], float32,
                                         (VariableDefinition.VARIABLE_NAME['time'], VariableDefinition.VARIABLE_NAME['latitude'], VariableDefinition.VARIABLE_NAME['longitude'],), fill_value=9.96921e+36)
        var.long_name = VariableDefinition.LONG_NAME['sea_surface_wave_peak_period']
        var.standard_name = VariableDefinition.STANDARD_NAME['sea_surface_wave_peak_period']
        var.units = VariableDefinition.CANONICAL_UNITS['sea_surface_wave_peak_period']

        logging.info('[DefaultWriter] Writing variable \'' + str(VariableDefinition.LONG_NAME['sea_surface_wave_peak_period']) + '\'')

        time_index = 0
        for time in self.coverage.read_axis_t():
            logging.debug(
                '[DefaultWriter] Writing variable \'' + str(VariableDefinition.LONG_NAME['sea_surface_wave_peak_period']) + '\' at time \'' + str(time) + '\'')
            var[time_index:time_index + 1, :] = self.coverage.read_variable_sea_surface_wave_peak_period_at_time(
                time)
            time_index += 1

    def write_variable_sea_surface_wave_from_direction(self):

        var = self.ncfile.createVariable(VariableDefinition.VARIABLE_NAME['sea_surface_wave_from_direction'], float32,
                                         (VariableDefinition.VARIABLE_NAME['time'], VariableDefinition.VARIABLE_NAME['latitude'], VariableDefinition.VARIABLE_NAME['longitude'],), fill_value=9.96921e+36)
        var.long_name = VariableDefinition.LONG_NAME['sea_surface_wave_from_direction']
        var.standard_name = VariableDefinition.STANDARD_NAME['sea_surface_wave_from_direction']
        var.units = VariableDefinition.CANONICAL_UNITS['sea_surface_wave_from_direction']

        logging.info('[DefaultWriter] Writing variable \'' + str(VariableDefinition.LONG_NAME['sea_surface_wave_from_direction']) + '\'')

        time_index = 0
        for time in self.coverage.read_axis_t():
            logging.debug(
                '[DefaultWriter] Writing variable \'' + str(VariableDefinition.LONG_NAME['sea_surface_wave_from_direction']) + '\' at time \'' + str(time) + '\'')
            var[time_index:time_index + 1, :] = self.coverage.read_variable_sea_surface_wave_from_direction_at_time(
                time)
            time_index += 1

    def write_variable_sea_surface_wave_to_direction(self):

        var = self.ncfile.createVariable(VariableDefinition.VARIABLE_NAME['sea_surface_wave_to_direction'], float32,
                                         (VariableDefinition.VARIABLE_NAME['time'], VariableDefinition.VARIABLE_NAME['latitude'], VariableDefinition.VARIABLE_NAME['longitude'],), fill_value=9.96921e+36)
        var.long_name = VariableDefinition.LONG_NAME['sea_surface_wave_to_direction']
        var.standard_name = VariableDefinition.STANDARD_NAME['sea_surface_wave_to_direction']
        var.units =VariableDefinition.CANONICAL_UNITS['sea_surface_wave_to_direction']

        logging.info('[DefaultWriter] Writing variable \'' + str(VariableDefinition.LONG_NAME['sea_surface_wave_to_direction']) + '\'')

        time_index = 0
        for time in self.coverage.read_axis_t():
            logging.debug(
                '[DefaultWriter] Writing variable \'' + str(VariableDefinition.LONG_NAME['sea_surface_wave_to_direction']) + '\' at time \'' + str(time) + '\'')
            var[time_index:time_index + 1, :] = self.coverage.read_variable_sea_surface_wave_to_direction_at_time(
                time)
            time_index += 1

    def write_variable_sea_surface_wave_stokes_drift_velocity(self):

        ucur = self.ncfile.createVariable(VariableDefinition.VARIABLE_NAME['eastward_sea_surface_wave_stokes_drift_velocity'], float32, (VariableDefinition.VARIABLE_NAME['time'],VariableDefinition.VARIABLE_NAME['latitude'], VariableDefinition.VARIABLE_NAME['longitude'],),fill_value=9.96921e+36)
        ucur.long_name = VariableDefinition.LONG_NAME['eastward_sea_surface_wave_stokes_drift_velocity']
        ucur.standard_name = VariableDefinition.STANDARD_NAME['eastward_sea_surface_wave_stokes_drift_velocity']
        ucur.units = VariableDefinition.CANONICAL_UNITS['eastward_sea_surface_wave_stokes_drift_velocity']

        vcur = self.ncfile.createVariable(VariableDefinition.VARIABLE_NAME['northward_sea_surface_wave_stokes_drift_velocity'], float32, (VariableDefinition.VARIABLE_NAME['time'], VariableDefinition.VARIABLE_NAME['latitude'], VariableDefinition.VARIABLE_NAME['longitude'],),fill_value=9.96921e+36)
        vcur.long_name = VariableDefinition.LONG_NAME['northward_sea_surface_wave_stokes_drift_velocity']
        vcur.standard_name = VariableDefinition.STANDARD_NAME['northward_sea_surface_wave_stokes_drift_velocity']
        vcur.units = VariableDefinition.CANONICAL_UNITS['northward_sea_surface_wave_stokes_drift_velocity']

        logging.info('[DefaultWriter] Writing variable \'Surface Stokes Drift Velocity\'\'')

        time_index=0
        for time in self.coverage.read_axis_t():

            logging.debug('[DefaultWriter] Writing variable \'Surface Stokes Drift Velocity\' at time \''+str(time)+'\'')

            cur = self.coverage.read_variable_sea_surface_wave_stokes_drift_velocity_at_time(time)

            ucur[time_index:time_index+1,:,:] = cur[0]
            vcur[time_index:time_index+1,:,:] = cur[1]
            time_index += 1

    def write_variable_atmosphere_momentum_flux_to_waves(self):

        ucur = self.ncfile.createVariable(VariableDefinition.VARIABLE_NAME['eastward_atmosphere_momentum_flux_to_waves'], float32, (VariableDefinition.VARIABLE_NAME['time'],VariableDefinition.VARIABLE_NAME['latitude'], VariableDefinition.VARIABLE_NAME['longitude'],),fill_value=9.96921e+36)
        ucur.long_name = VariableDefinition.LONG_NAME['eastward_atmosphere_momentum_flux_to_waves']
        ucur.standard_name = VariableDefinition.STANDARD_NAME['eastward_atmosphere_momentum_flux_to_waves']
        ucur.units = VariableDefinition.CANONICAL_UNITS['eastward_atmosphere_momentum_flux_to_waves']

        vcur = self.ncfile.createVariable(VariableDefinition.VARIABLE_NAME['northward_atmosphere_momentum_flux_to_waves'], float32, (VariableDefinition.VARIABLE_NAME['time'], VariableDefinition.VARIABLE_NAME['latitude'], VariableDefinition.VARIABLE_NAME['longitude'],),fill_value=9.96921e+36)
        vcur.long_name = VariableDefinition.LONG_NAME['northward_atmosphere_momentum_flux_to_waves']
        vcur.standard_name = VariableDefinition.STANDARD_NAME['northward_atmosphere_momentum_flux_to_waves']
        vcur.units =VariableDefinition.CANONICAL_UNITS['northward_atmosphere_momentum_flux_to_waves']

        logging.info('[DefaultWriter] Writing variable \'Atmosphere Momentum Flux to Waves\'\'')

        time_index=0
        for time in self.coverage.read_axis_t():

            logging.debug('[DefaultWriter] Writing variable \'Atmosphere Momentum Flux to Waves\' at time \''+str(time)+'\'')

            cur = self.coverage.read_variable_atmosphere_momentum_flux_to_waves_at_time(time)

            ucur[time_index:time_index+1,:,:] = cur[0]
            vcur[time_index:time_index+1,:,:] = cur[1]
            time_index += 1

    def write_variable_waves_momentum_flux_to_ocean(self):

        ucur = self.ncfile.createVariable(VariableDefinition.VARIABLE_NAME['eastward_waves_momentum_flux_to_ocean'], float32, (VariableDefinition.VARIABLE_NAME['time'],VariableDefinition.VARIABLE_NAME['latitude'], VariableDefinition.VARIABLE_NAME['longitude'],),fill_value=9.96921e+36)
        ucur.long_name = VariableDefinition.LONG_NAME['eastward_waves_momentum_flux_to_ocean']
        ucur.standard_name = VariableDefinition.STANDARD_NAME['eastward_waves_momentum_flux_to_ocean']
        ucur.units = VariableDefinition.CANONICAL_UNITS['eastward_waves_momentum_flux_to_ocean']

        vcur = self.ncfile.createVariable(VariableDefinition.VARIABLE_NAME['northward_waves_momentum_flux_to_ocean'], float32, (VariableDefinition.VARIABLE_NAME['time'], VariableDefinition.VARIABLE_NAME['latitude'], VariableDefinition.VARIABLE_NAME['longitude'],),fill_value=9.96921e+36)
        vcur.long_name = VariableDefinition.LONG_NAME['northward_waves_momentum_flux_to_ocean']
        vcur.standard_name = VariableDefinition.STANDARD_NAME['northward_waves_momentum_flux_to_ocean']
        vcur.units = VariableDefinition.CANONICAL_UNITS['northward_waves_momentum_flux_to_ocean']

        logging.info('[DefaultWriter] Writing variable \'Waves Momentum Flux To Ocean\'\'')

        time_index=0
        for time in self.coverage.read_axis_t():

            logging.debug('[DefaultWriter] Writing variable \'Waves Momentum Flux To Ocean\' at time \''+str(time)+'\'')

            cur = self.coverage.read_variable_waves_momentum_flux_to_ocean_at_time(time)

            ucur[time_index:time_index+1,:,:] = cur[0]
            vcur[time_index:time_index+1,:,:] = cur[1]
            time_index += 1

    def write_variable_radiation_pressure_bernouilli_head(self):

        var = self.ncfile.createVariable(VariableDefinition.VARIABLE_NAME['radiation_pressure_bernouilli_head'], float32,
                                         (VariableDefinition.VARIABLE_NAME['time'], VariableDefinition.VARIABLE_NAME['latitude'], VariableDefinition.VARIABLE_NAME['longitude'],), fill_value=9.96921e+36)
        var.long_name = VariableDefinition.LONG_NAME['radiation_pressure_bernouilli_head']
        var.standard_name = VariableDefinition.STANDARD_NAME['radiation_pressure_bernouilli_head']
        var.units = VariableDefinition.CANONICAL_UNITS['radiation_pressure_bernouilli_head']

        logging.info('[DefaultWriter] Writing variable \'' + str(VariableDefinition.LONG_NAME['radiation_pressure_bernouilli_head']) + '\'')

        time_index = 0
        for time in self.coverage.read_axis_t():
            logging.debug(
                '[DefaultWriter] Writing variable \'' + str(VariableDefinition.LONG_NAME['radiation_pressure_bernouilli_head']) + '\' at time \'' + str(time) + '\'')
            var[time_index:time_index + 1,
            :] = self.coverage.read_variable_radiation_pressure_bernouilli_head_at_time(time)
            time_index += 1

    def write_variable_sea_surface_wave_energy_dissipation_at_ground_level(self):

        var = self.ncfile.createVariable(VariableDefinition.VARIABLE_NAME['sea_surface_wave_energy_dissipation_at_ground_level'], float32,
                                         (VariableDefinition.VARIABLE_NAME['time'], VariableDefinition.VARIABLE_NAME['latitude'], VariableDefinition.VARIABLE_NAME['longitude'],), fill_value=9.96921e+36)
        var.long_name = VariableDefinition.LONG_NAME['sea_surface_wave_energy_dissipation_at_ground_level']
        var.standard_name = VariableDefinition.STANDARD_NAME['sea_surface_wave_energy_dissipation_at_ground_level']
        var.units = VariableDefinition.CANONICAL_UNITS['sea_surface_wave_energy_dissipation_at_ground_level']

        logging.info('[DefaultWriter] Writing variable \'' + str(VariableDefinition.LONG_NAME['sea_surface_wave_energy_dissipation_at_ground_level']) + '\'')

        time_index = 0
        for time in self.coverage.read_axis_t():
            logging.debug(
                '[DefaultWriter] Writing variable \'' + str(VariableDefinition.LONG_NAME['sea_surface_wave_energy_dissipation_at_ground_level']) + '\' at time \'' + str(time) + '\'')
            var[time_index:time_index + 1,
            :] = self.coverage.read_variable_sea_surface_wave_energy_dissipation_at_ground_level_at_time(time)
            time_index += 1

    def write_variable_sea_surface_wave_energy_flux_to_ocean(self):

        var = self.ncfile.createVariable(VariableDefinition.VARIABLE_NAME['sea_surface_wave_energy_flux_to_ocean'], float32,
                                         (VariableDefinition.VARIABLE_NAME['time'], VariableDefinition.VARIABLE_NAME['latitude'], VariableDefinition.VARIABLE_NAME['longitude'],), fill_value=9.96921e+36)
        var.long_name = VariableDefinition.LONG_NAME['sea_surface_wave_energy_flux_to_ocean']
        var.standard_name = VariableDefinition.STANDARD_NAME['sea_surface_wave_energy_flux_to_ocean']
        var.units = VariableDefinition.CANONICAL_UNITS['sea_surface_wave_energy_flux_to_ocean']

        logging.info('[DefaultWriter] Writing variable \'' + str(VariableDefinition.LONG_NAME['sea_surface_wave_energy_flux_to_ocean']) + '\'')

        time_index = 0
        for time in self.coverage.read_axis_t():
            logging.debug(
                '[DefaultWriter] Writing variable \'' + str(VariableDefinition.LONG_NAME['sea_surface_wave_energy_flux_to_ocean']) + '\' at time \'' + str(time) + '\'')
            var[time_index:time_index + 1,
            :] = self.coverage.read_variable_sea_surface_wave_energy_flux_to_ocean_at_time(time)
            time_index += 1

    # METEO
    def write_variable_wind_10m(self):

        ucur = self.ncfile.createVariable(VariableDefinition.VARIABLE_NAME['eastward_wind_10m'], float32, (VariableDefinition.VARIABLE_NAME['time'],VariableDefinition.VARIABLE_NAME['latitude'], VariableDefinition.VARIABLE_NAME['longitude'],),fill_value=9.96921e+36)
        ucur.long_name = VariableDefinition.LONG_NAME['eastward_wind_10m']
        ucur.standard_name = VariableDefinition.STANDARD_NAME['eastward_wind_10m']
        ucur.units = VariableDefinition.CANONICAL_UNITS['eastward_wind_10m']

        vcur = self.ncfile.createVariable(VariableDefinition.VARIABLE_NAME['northward_wind_10m'], float32, (VariableDefinition.VARIABLE_NAME['time'], VariableDefinition.VARIABLE_NAME['latitude'], VariableDefinition.VARIABLE_NAME['longitude'],),fill_value=9.96921e+36)
        vcur.long_name = VariableDefinition.LONG_NAME['northward_wind_10m']
        vcur.standard_name = VariableDefinition.STANDARD_NAME['northward_wind_10m']
        vcur.units = VariableDefinition.CANONICAL_UNITS['northward_wind_10m']

        logging.info('[DefaultWriter] Writing variable \'Wind 10m\'\'')

        time_index=0
        for time in self.coverage.read_axis_t():

            logging.debug('[DefaultWriter] Writing variable \'Wind 10m\' at time \''+str(time)+'\'')

            cur = self.coverage.read_variable_wind_10m_at_time(time)

            ucur[time_index:time_index+1,:,:] = cur[0]
            vcur[time_index:time_index+1,:,:] = cur[1]
            time_index += 1

    def write_variable_wind_stress(self):

        ucur = self.ncfile.createVariable(VariableDefinition.VARIABLE_NAME['eastward_wind_stress'], float32, (VariableDefinition.VARIABLE_NAME['time'], VariableDefinition.VARIABLE_NAME['latitude'], VariableDefinition.VARIABLE_NAME['longitude'],),
                                          fill_value=9.96921e+36)
        ucur.long_name = VariableDefinition.LONG_NAME['eastward_wind_stress']
        ucur.standard_name = VariableDefinition.STANDARD_NAME['eastward_wind_stress']
        ucur.units = VariableDefinition.CANONICAL_UNITS['eastward_wind_stress']

        vcur = self.ncfile.createVariable(VariableDefinition.VARIABLE_NAME['northward_wind_stress'], float32, (VariableDefinition.VARIABLE_NAME['time'], VariableDefinition.VARIABLE_NAME['latitude'], VariableDefinition.VARIABLE_NAME['longitude'],),
                                          fill_value=9.96921e+36)
        vcur.long_name = VariableDefinition.LONG_NAME['northward_wind_stress']
        vcur.standard_name = VariableDefinition.STANDARD_NAME['northward_wind_stress']
        vcur.units =VariableDefinition.CANONICAL_UNITS['northward_wind_stress']

        logging.info('[DefaultWriter] Writing variable \'Wind Stress\'\'')

        time_index = 0
        for time in self.coverage.read_axis_t():
            logging.debug('[DefaultWriter] Writing variable \'Wind Stress\' at time \'' + str(time) + '\'')

            cur = self.coverage.read_variable_wind_stress_at_time(time)

            ucur[time_index:time_index + 1, :, :] = cur[0]
            vcur[time_index:time_index + 1, :, :] = cur[1]
            time_index += 1



