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
from utils.VariableUnits import VariableUnits
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
        self.ncfile.createDimension('latitude', self.coverage.get_y_size())
        self.ncfile.createDimension('longitude', self.coverage.get_x_size())

        if self.coverage.is_regular_grid()==True:

            # variables
            latitudes = self.ncfile.createVariable('latitude', float32, ('latitude',))
            latitudes.long_name = "latitude" ;
            latitudes.standard_name = "latitude" ;
            latitudes.valid_min = -90.;
            latitudes.valid_max = 90. ;
            latitudes.axis = "Y" ;
            latitudes.units = VariableUnits.CANONICAL_UNITS[latitudes.standard_name];

            longitudes = self.ncfile.createVariable('longitude', float32, ('longitude',))
            longitudes.long_name = "longitude" ;
            longitudes.standard_name = "longitude" ;
            longitudes.valid_min = -180. ;
            longitudes.valid_max = 180. ;
            longitudes.axis = "X" ;
            longitudes.units = VariableUnits.CANONICAL_UNITS[longitudes.standard_name];

            # data
            latitudes[:] = self.coverage.read_axis_y();
            longitudes[:] = self.coverage.read_axis_x();

        else:

            latitudes = self.ncfile.createVariable('latitude', float32, ('latitude','longitude',))
            latitudes.long_name = "latitude" ;
            latitudes.standard_name = "latitude" ;
            latitudes.valid_min = -90.;
            latitudes.valid_max = 90. ;
            latitudes.axis = "Y" ;
            latitudes.units = VariableUnits.CANONICAL_UNITS[latitudes.standard_name];

            longitudes = self.ncfile.createVariable('longitude', float32, ('latitude','longitude',))
            longitudes.long_name = "longitude" ;
            longitudes.standard_name = "longitude" ;
            longitudes.valid_min = -180. ;
            longitudes.valid_max = 180. ;
            longitudes.axis = "X" ;
            longitudes.units = VariableUnits.CANONICAL_UNITS[longitudes.standard_name];

            # data
            latitudes[:,:] = self.coverage.read_axis_y();
            longitudes[:,:] = self.coverage.read_axis_x();


        if(isinstance(self.coverage, TimeCoverage) or isinstance(self.coverage, TimeLevelCoverage)):

            self.ncfile.createDimension('time', None)
            times = self.ncfile.createVariable('time', float64, ('time',))
            times.units= 'seconds since 1970-01-01 00:00:00'
            times.calendar= 'gregorian'
            times.standard_name= 'time'
            times.axis='T'
            times.conventions = "UTC time"

            times[:] = date2num(self.coverage.read_axis_t(), units = times.units, calendar = times.calendar)

        if(isinstance(self.coverage, LevelCoverage) or isinstance(self.coverage, TimeLevelCoverage)):

            if self.coverage.is_sigma_coordinate():
                raise ValueError("This writer supports only Coverage with a regular vertical axis.")

            self.ncfile.createDimension('depth', self.coverage.get_z_size())
            levels = self.ncfile.createVariable('depth', float64, ('depth',))
            levels.standard_name= 'depth'
            levels.long_name="Positive depth"
            levels.axis='Z'
            levels.units = VariableUnits.CANONICAL_UNITS[levels.standard_name];

            levels[:] = self.coverage.read_axis_z()

        
    def close(self):
        self.ncfile.close()

    # HYDRO
    def write_variable_mesh_size(self):

        var = self.ncfile.createVariable('mesh_size', float32, ('latitude', 'longitude',), fill_value=9.96921e+36)
        var.long_name = "Mesh size";
        var.standard_name = "mesh_size";
        var.units = VariableUnits.CANONICAL_UNITS[var.standard_name];

        logging.info('[DefaultWriter] Writing variable \'' + str(var.long_name) + '\'')

        var[:] = self.coverage.read_variable_mesh_size()

    def write_variable_2D_sea_binary_mask(self):

        var = self.ncfile.createVariable('2d_sea_binary_mask', int, ('latitude', 'longitude',), fill_value=-9999)
        var.long_name = "Sea Binary Mask";
        var.standard_name = "sea_binary_mask";
        var.units = VariableUnits.CANONICAL_UNITS[var.standard_name];
        var.comment = "1 = sea, 0 = land"

        logging.info('[DefaultWriter] Writing variable \'' + str(var.long_name) + '\'')
        var[:] = self.coverage.read_variable_2D_sea_binary_mask()

    def write_variable_2D_wet_binary_mask(self):

        var = self.ncfile.createVariable('wet_binary_mask', float32, ('time', 'latitude', 'longitude',),
                                         fill_value=9.96921e+36)
        var.long_name = "wet_binary_mask";
        var.standard_name = "wet_binary_mask";
        var.units = VariableUnits.CANONICAL_UNITS[var.standard_name];

        logging.info('[DefaultWriter] Writing variable \'' + str(var.long_name) + '\'')

        time_index = 0
        for time in self.coverage.read_axis_t():
            logging.debug(
                '[DefaultWriter] Writing variable \'' + str(var.long_name) + '\' at time \'' + str(time) + '\'')
            var[time_index:time_index + 1, :] = self.coverage.read_variable_2D_wet_binary_mask_at_time(time)
            time_index += 1

    def write_variable_bathymetry(self):

        var = self.ncfile.createVariable('bathymetry', float32, ('latitude', 'longitude',),fill_value=9.96921e+36)
        var.long_name = "Bathymetry" ;
        var.standard_name = "bathymetry" ;
        var.units = VariableUnits.CANONICAL_UNITS[var.standard_name];

        logging.info('[DefaultWriter] Writing variable \'' + str(var.long_name) + '\'')

        var[:] = self.coverage.read_variable_bathymetry()
            
    def write_variable_sea_surface_height_above_mean_sea_level(self):
            
        var = self.ncfile.createVariable('sea_surface_height_above_mean_sea_level', float32, ('time', 'latitude', 'longitude',),fill_value=9.96921e+36)
        var.long_name = "Sea Surface Height Above Mean Sea Level" ;
        var.standard_name = "sea_surface_height_above_mean_sea_level" ;
        var.units = VariableUnits.CANONICAL_UNITS[var.standard_name];

        logging.info('[DefaultWriter] Writing variable \''+str(var.long_name)+'\'')
        
        time_index=0
        for time in self.coverage.read_axis_t():
            logging.debug('[DefaultWriter] Writing variable \''+str(var.long_name)+'\' at time \''+str(time)+'\'')
            var[time_index:time_index+1,:] = self.coverage.read_variable_sea_surface_height_above_mean_sea_level_at_time(time)
            time_index += 1

    def write_variable_sea_surface_height_above_goeid(self):

        var = self.ncfile.createVariable('sea_surface_height_above_goeid', float32,
                                         ('time', 'latitude', 'longitude',), fill_value=9.96921e+36)
        var.long_name = "Sea Surface Height Above Geoid";
        var.standard_name = "sea_surface_height_above_goeid";
        var.units = VariableUnits.CANONICAL_UNITS[var.standard_name];

        logging.info('[DefaultWriter] Writing variable \'' + str(var.long_name) + '\'')

        time_index = 0
        for time in self.coverage.read_axis_t():
            logging.debug(
                '[DefaultWriter] Writing variable \'' + str(var.long_name) + '\' at time \'' + str(time) + '\'')
            var[time_index:time_index + 1,
            :] = self.coverage.read_variable_sea_surface_height_above_goeid_at_time(time)
            time_index += 1

    def write_variable_baroclinic_sea_water_velocity(self):

        ucur = self.ncfile.createVariable('baroclinic_eastward_sea_water_velocity', float32,
                                          ('time', 'depth', 'latitude', 'longitude',), fill_value=9.96921e+36)
        ucur.long_name = "Baroclinic Eastward Sea Water Velocity";
        ucur.standard_name = "baroclinic_eastward_sea_water_velocity";
        ucur.units = VariableUnits.CANONICAL_UNITS[ucur.standard_name];
        ucur.comment = "cur=sqrt(U**2+V**2)";

        vcur = self.ncfile.createVariable('baroclinic_northward_sea_water_velocity', float32,
                                          ('time', 'depth', 'latitude', 'longitude',), fill_value=9.96921e+36)
        vcur.long_name = "Baroclinic Northward Sea Water Velocity";
        vcur.standard_name = "baroclinic_northward_sea_water_velocity";
        vcur.units = VariableUnits.CANONICAL_UNITS[vcur.standard_name];
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
        ucur = self.ncfile.createVariable('barotropic_eastward_sea_water_velocity', float32,
                                          ('time', 'latitude', 'longitude',), fill_value=9.96921e+36)
        ucur.long_name = "Barotropic Eastward Sea Water Velocity";
        ucur.standard_name = "barotropic_eastward_sea_water_velocity";
        ucur.units = VariableUnits.CANONICAL_UNITS[ucur.standard_name];
        ucur.comment = "cur=sqrt(U**2+V**2)";

        vcur = self.ncfile.createVariable('barotropic_northward_sea_water_velocity', float32,
                                          ('time', 'latitude', 'longitude',), fill_value=9.96921e+36)
        vcur.long_name = "Barotropic Northward Sea Water Velocity";
        vcur.standard_name = "barotropic_northward_sea_water_velocity";
        vcur.units = VariableUnits.CANONICAL_UNITS[vcur.standard_name];
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

        var = self.ncfile.createVariable('sea_surface_temperature', float32, ('time', 'latitude', 'longitude',),
                                         fill_value=9.96921e+36)
        var.long_name = "Sea Surface Temperature";
        var.standard_name = "sea_surface_temperature";
        var.units = VariableUnits.CANONICAL_UNITS[var.standard_name];

        logging.info('[DefaultWriter] Writing variable \'' + str(var.long_name) + '\'')

        time_index = 0
        for time in self.coverage.read_axis_t():
            logging.debug(
                '[DefaultWriter] Writing variable \'' + str(var.long_name) + '\' at time \'' + str(time) + '\'')
            var[time_index:time_index + 1, :] = self.coverage.read_variable_sea_surface_temperature_at_time(time)
            time_index += 1

    def write_variable_sea_water_temperature(self):

        var = self.ncfile.createVariable('sea_water_temperature', float32,
                                         ('time', 'depth', 'latitude', 'longitude',), fill_value=9.96921e+36)
        var.long_name = "Sea Water Temperature";
        var.standard_name = "sea_water_temperature";
        var.units = VariableUnits.CANONICAL_UNITS[var.standard_name];

        logging.info('[DefaultWriter] Writing variable \'' + str(var.long_name) + '\'')

        time_index = 0
        for time in self.coverage.read_axis_t():

            logging.debug(
                '[DefaultWriter] Writing variable \'' + str(var.long_name) + '\' at time \'' + str(time) + '\'')

            level_index = 0
            for level in self.coverage.read_axis_z():
                var[time_index:time_index + 1, level_index:level_index + 1, :,
                :] = self.coverage.read_variable_sea_water_temperature_at_time_and_depth(time, level)
                level_index += 1

            time_index += 1

    def write_variable_sea_water_salinity(self):

        var = self.ncfile.createVariable('sea_water_salinity', float32,
                                         ('time', 'depth', 'latitude', 'longitude',), fill_value=9.96921e+36)
        var.long_name = "Sea Water Salinity";
        var.standard_name = "sea_water_salinity";
        var.units = VariableUnits.CANONICAL_UNITS[var.standard_name];

        logging.info('[DefaultWriter] Writing variable \'' + str(var.long_name) + '\'')

        time_index = 0
        for time in self.coverage.read_axis_t():

            logging.debug(
                '[DefaultWriter] Writing variable \'' + str(var.long_name) + '\' at time \'' + str(time) + '\'')

            level_index = 0
            for level in self.coverage.read_axis_z():
                var[time_index:time_index + 1, level_index:level_index + 1, :,
                :] = self.coverage.read_variable_sea_water_salinity_at_time_and_depth(time, level)
                level_index += 1

            time_index += 1

    # WAVES
    def write_variable_sea_surface_wave_significant_height(self):

        var = self.ncfile.createVariable('sea_surface_wave_significant_height', float32, ('time', 'latitude', 'longitude',),fill_value=9.96921e+36)
        var.long_name = "Sea Surface Wave Significant Height" ;
        var.standard_name = "sea_surface_wave_significant_height" ;
        var.units = VariableUnits.CANONICAL_UNITS[var.standard_name];

        logging.info('[DefaultWriter] Writing variable \'' + str(var.long_name) + '\'')

        time_index=0
        for time in self.coverage.read_axis_t():
            logging.debug(
                '[DefaultWriter] Writing variable \'' + str(var.long_name) + '\' at time \'' + str(time) + '\'')
            var[time_index:time_index+1,:] = self.coverage.read_variable_sea_surface_wave_significant_height_at_time(time)
            time_index += 1

    def write_variable_sea_surface_wave_breaking_height(self):

        var = self.ncfile.createVariable('sea_surface_wave_breaking_height', float32,
                                         ('time', 'latitude', 'longitude',), fill_value=9.96921e+36)
        var.long_name = "Sea Surface Wave Breaking Height";
        var.standard_name = "sea_surface_wave_breaking_height";
        var.units = VariableUnits.CANONICAL_UNITS[var.standard_name];

        logging.info('[DefaultWriter] Writing variable \'' + str(var.long_name) + '\'')

        time_index = 0
        for time in self.coverage.read_axis_t():
            logging.debug(
                '[DefaultWriter] Writing variable \'' + str(var.long_name) + '\' at time \'' + str(time) + '\'')
            var[time_index:time_index + 1,
            :] = self.coverage.read_variable_sea_surface_wave_breaking_height_at_time(time)
            time_index += 1

    def write_variable_sea_surface_wave_mean_period(self):

        var = self.ncfile.createVariable('sea_surface_wave_mean_period', float32,
                                         ('time', 'latitude', 'longitude',), fill_value=9.96921e+36)
        var.long_name = "Sea Surface Wave Mean_Period";
        var.standard_name = "sea_surface_wave_mean_period";
        var.units = VariableUnits.CANONICAL_UNITS[var.standard_name];

        logging.info('[DefaultWriter] Writing variable \'' + str(var.long_name) + '\'')

        time_index = 0
        for time in self.coverage.read_axis_t():
            logging.debug(
                '[DefaultWriter] Writing variable \'' + str(var.long_name) + '\' at time \'' + str(time) + '\'')
            var[time_index:time_index + 1, :] = self.coverage.read_variable_sea_surface_wave_mean_period_at_time(
                time)
            time_index += 1

    def write_variable_sea_surface_wave_peak_period(self):

        var = self.ncfile.createVariable('sea_surface_wave_peak_period', float32,
                                         ('time', 'latitude', 'longitude',), fill_value=9.96921e+36)
        var.long_name = "Sea Surface Wave Peak_Period";
        var.standard_name = "sea_surface_wave_peak_period";
        var.units = VariableUnits.CANONICAL_UNITS[var.standard_name];

        logging.info('[DefaultWriter] Writing variable \'' + str(var.long_name) + '\'')

        time_index = 0
        for time in self.coverage.read_axis_t():
            logging.debug(
                '[DefaultWriter] Writing variable \'' + str(var.long_name) + '\' at time \'' + str(time) + '\'')
            var[time_index:time_index + 1, :] = self.coverage.read_variable_sea_surface_wave_peak_period_at_time(
                time)
            time_index += 1

    def write_variable_sea_surface_wave_from_direction(self):

        var = self.ncfile.createVariable('sea_surface_wave_from_direction', float32,
                                         ('time', 'latitude', 'longitude',), fill_value=9.96921e+36)
        var.long_name = "Sea Surface From Direction";
        var.standard_name = "sea_surface_wave_from_direction";
        var.units = VariableUnits.CANONICAL_UNITS[var.standard_name];

        logging.info('[DefaultWriter] Writing variable \'' + str(var.long_name) + '\'')

        time_index = 0
        for time in self.coverage.read_axis_t():
            logging.debug(
                '[DefaultWriter] Writing variable \'' + str(var.long_name) + '\' at time \'' + str(time) + '\'')
            var[time_index:time_index + 1, :] = self.coverage.read_variable_sea_surface_wave_from_direction_at_time(
                time)
            time_index += 1

    def write_variable_sea_surface_wave_to_direction(self):

        var = self.ncfile.createVariable('sea_surface_wave_to_direction', float32,
                                         ('time', 'latitude', 'longitude',), fill_value=9.96921e+36)
        var.long_name = "Sea Surface To Direction";
        var.standard_name = "sea_surface_wave_to_direction";
        var.units = VariableUnits.CANONICAL_UNITS[var.standard_name];

        logging.info('[DefaultWriter] Writing variable \'' + str(var.long_name) + '\'')

        time_index = 0
        for time in self.coverage.read_axis_t():
            logging.debug(
                '[DefaultWriter] Writing variable \'' + str(var.long_name) + '\' at time \'' + str(time) + '\'')
            var[time_index:time_index + 1, :] = self.coverage.read_variable_sea_surface_wave_to_direction_at_time(
                time)
            time_index += 1

    def write_variable_sea_surface_wave_stokes_drift_velocity(self):

        ucur = self.ncfile.createVariable('eastward_sea_surface_wave_stokes_drift_velocity', float32, ('time','latitude', 'longitude',),fill_value=9.96921e+36)
        ucur.long_name = "Eastward Sea Surface Wave Stokes Dift Velocity" ;
        ucur.standard_name = "eastward_sea_surface_wave_stokes_drift_velocity" ;
        ucur.units = VariableUnits.CANONICAL_UNITS[ucur.standard_name];

        vcur = self.ncfile.createVariable('northward_sea_surface_wave_stokes_drift_velocity', float32, ('time', 'latitude', 'longitude',),fill_value=9.96921e+36)
        vcur.long_name = "Northward Sea Surface Wave Stokes Drift Velocity" ;
        vcur.standard_name = "northward_sea_surface_wave_stokes_drift_velocity" ;
        vcur.units = VariableUnits.CANONICAL_UNITS[vcur.standard_name];

        logging.info('[DefaultWriter] Writing variable \'Surface Stokes Drift Velocity\'\'')

        time_index=0
        for time in self.coverage.read_axis_t():

            logging.debug('[DefaultWriter] Writing variable \'Surface Stokes Drift Velocity\' at time \''+str(time)+'\'')

            cur = self.coverage.read_variable_sea_surface_wave_stokes_drift_velocity_at_time(time)

            ucur[time_index:time_index+1,:,:] = cur[0]
            vcur[time_index:time_index+1,:,:] = cur[1]
            time_index += 1

    def write_variable_atmosphere_momentum_flux_to_waves(self):

        ucur = self.ncfile.createVariable('eastward_atmosphere_momentum_flux_to_waves', float32, ('time','latitude', 'longitude',),fill_value=9.96921e+36)
        ucur.long_name = "Eastward Atmosphere Momentum Flux to Waves" ;
        ucur.standard_name = "eastward_atmosphere_momentum_flux_to_waves" ;
        ucur.units = VariableUnits.CANONICAL_UNITS[ucur.standard_name];

        vcur = self.ncfile.createVariable('northward_atmosphere_momentum_flux_to_waves', float32, ('time', 'latitude', 'longitude',),fill_value=9.96921e+36)
        vcur.long_name = "Northward Atmosphere Momentum Flux to Waves" ;
        vcur.standard_name = "northward_atmosphere_momentum_flux_to_waves" ;
        vcur.units = VariableUnits.CANONICAL_UNITS[vcur.standard_name];

        logging.info('[DefaultWriter] Writing variable \'Atmosphere Momentum Flux to Waves\'\'')

        time_index=0
        for time in self.coverage.read_axis_t():

            logging.debug('[DefaultWriter] Writing variable \'Atmosphere Momentum Flux to Waves\' at time \''+str(time)+'\'')

            cur = self.coverage.read_variable_atmosphere_momentum_flux_to_waves_at_time(time)

            ucur[time_index:time_index+1,:,:] = cur[0]
            vcur[time_index:time_index+1,:,:] = cur[1]
            time_index += 1

    def write_variable_waves_momentum_flux_to_ocean(self):

        ucur = self.ncfile.createVariable('eastward_waves_momentum_flux_to_ocean', float32, ('time','latitude', 'longitude',),fill_value=9.96921e+36)
        ucur.long_name = "Eastward Waves Momentum Flux To Ocean" ;
        ucur.standard_name = "eastward_waves_momentum_flux_to_ocean" ;
        ucur.units = VariableUnits.CANONICAL_UNITS[ucur.standard_name];


        vcur = self.ncfile.createVariable('northward_waves_momentum_flux_to_ocean', float32, ('time', 'latitude', 'longitude',),fill_value=9.96921e+36)
        vcur.long_name = "Northward Waves Momentum Flux To Ocean" ;
        vcur.standard_name = "northward_waves_momentum_flux_to_ocean" ;
        vcur.units = VariableUnits.CANONICAL_UNITS[vcur.standard_name];

        logging.info('[DefaultWriter] Writing variable \'Waves Momentum Flux To Ocean\'\'')

        time_index=0
        for time in self.coverage.read_axis_t():

            logging.debug('[DefaultWriter] Writing variable \'Waves Momentum Flux To Ocean\' at time \''+str(time)+'\'')

            cur = self.coverage.read_variable_waves_momentum_flux_to_ocean_at_time(time)

            ucur[time_index:time_index+1,:,:] = cur[0]
            vcur[time_index:time_index+1,:,:] = cur[1]
            time_index += 1

    def write_variable_radiation_pressure_bernouilli_head(self):

        var = self.ncfile.createVariable('radiation_pressure_bernouilli_head', float32,
                                         ('time', 'latitude', 'longitude',), fill_value=9.96921e+36)
        var.long_name = "Radiation Pressure (Bernouilli Head)";
        var.standard_name = "radiation_pressure_bernouilli_head";
        var.units = VariableUnits.CANONICAL_UNITS[var.standard_name];

        logging.info('[DefaultWriter] Writing variable \'' + str(var.long_name) + '\'')

        time_index = 0
        for time in self.coverage.read_axis_t():
            logging.debug(
                '[DefaultWriter] Writing variable \'' + str(var.long_name) + '\' at time \'' + str(time) + '\'')
            var[time_index:time_index + 1,
            :] = self.coverage.read_variable_radiation_pressure_bernouilli_head_at_time(time)
            time_index += 1

    def write_variable_sea_surface_wave_energy_dissipation_at_ground_level(self):

        var = self.ncfile.createVariable('sea_surface_wave_energy_dissipation_at_ground_level', float32,
                                         ('time', 'latitude', 'longitude',), fill_value=9.96921e+36)
        var.long_name = "Sea Surface Wave Energy Dissipation At Ground Level";
        var.standard_name = "sea_surface_wave_energy_dissipation_at_ground_level";
        var.units = VariableUnits.CANONICAL_UNITS[var.standard_name];

        logging.info('[DefaultWriter] Writing variable \'' + str(var.long_name) + '\'')

        time_index = 0
        for time in self.coverage.read_axis_t():
            logging.debug(
                '[DefaultWriter] Writing variable \'' + str(var.long_name) + '\' at time \'' + str(time) + '\'')
            var[time_index:time_index + 1,
            :] = self.coverage.read_variable_sea_surface_wave_energy_dissipation_at_ground_level_at_time(time)
            time_index += 1

    def write_variable_sea_surface_wave_energy_flux_to_ocean(self):

        var = self.ncfile.createVariable('sea_surface_wave_energy_flux_to_ocean', float32,
                                         ('time', 'latitude', 'longitude',), fill_value=9.96921e+36)
        var.long_name = "Sea Surface Wave Energy Flux To Ocean";
        var.standard_name = "sea_surface_wave_energy_flux_to_ocean";
        var.units = VariableUnits.CANONICAL_UNITS[var.standard_name];

        logging.info('[DefaultWriter] Writing variable \'' + str(var.long_name) + '\'')

        time_index = 0
        for time in self.coverage.read_axis_t():
            logging.debug(
                '[DefaultWriter] Writing variable \'' + str(var.long_name) + '\' at time \'' + str(time) + '\'')
            var[time_index:time_index + 1,
            :] = self.coverage.read_variable_sea_surface_wave_energy_flux_to_ocean_at_time(time)
            time_index += 1

    # METEO
    def write_variable_wind_10m(self):

        ucur = self.ncfile.createVariable('eastward_wind_10m', float32, ('time','latitude', 'longitude',),fill_value=9.96921e+36)
        ucur.long_name = "Eastward Wind 10m" ;
        ucur.standard_name = "eastward_wind_10m" ;
        ucur.units = VariableUnits.CANONICAL_UNITS[ucur.standard_name];

        vcur = self.ncfile.createVariable('northward_wind_10m', float32, ('time', 'latitude', 'longitude',),fill_value=9.96921e+36)
        vcur.long_name = "Northward Wind 10m" ;
        vcur.standard_name = "northward_wind_10m" ;
        vcur.units = VariableUnits.CANONICAL_UNITS[vcur.standard_name];

        logging.info('[DefaultWriter] Writing variable \'Wind 10m\'\'')

        time_index=0
        for time in self.coverage.read_axis_t():

            logging.debug('[DefaultWriter] Writing variable \'Wind 10m\' at time \''+str(time)+'\'')

            cur = self.coverage.read_variable_wind_10m_at_time(time)

            ucur[time_index:time_index+1,:,:] = cur[0]
            vcur[time_index:time_index+1,:,:] = cur[1]
            time_index += 1

    def write_variable_wind_stress(self):

        ucur = self.ncfile.createVariable('eastward_wind_stress', float32, ('time', 'latitude', 'longitude',),
                                          fill_value=9.96921e+36)
        ucur.long_name = "Eastward Wind Stress";
        ucur.standard_name = "eastward_wind_stress";
        ucur.units = VariableUnits.CANONICAL_UNITS[ucur.standard_name];

        vcur = self.ncfile.createVariable('northward_wind_stress', float32, ('time', 'latitude', 'longitude',),
                                          fill_value=9.96921e+36)
        vcur.long_name = "Northward Wind Stress";
        vcur.standard_name = "northward_wind_stress";
        vcur.units = VariableUnits.CANONICAL_UNITS[vcur.standard_name];

        logging.info('[DefaultWriter] Writing variable \'Wind Stress\'\'')

        time_index = 0
        for time in self.coverage.read_axis_t():
            logging.debug('[DefaultWriter] Writing variable \'Wind Stress\' at time \'' + str(time) + '\'')

            cur = self.coverage.read_variable_wind_stress_at_time(time)

            ucur[time_index:time_index + 1, :, :] = cur[0]
            vcur[time_index:time_index + 1, :, :] = cur[1]
            time_index += 1



