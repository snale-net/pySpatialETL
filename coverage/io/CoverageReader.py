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
import os

class CoverageReader(object):
    
    def __init__(self, myFile):

        if myFile is not None and not os.path.exists(myFile):
            raise IOError(str(myFile) + " doesn't exists. Abort")

        self.filename = myFile;

    def read_axis_x(self):
        raise NotImplementedError(str(type(self))+" don't have implemented the function 'read_axis_x()'.")

    def read_axis_y(self):
        raise NotImplementedError(str(type(self)) + " don't have implemented the function 'read_axis_y()'.")

    def read_axis_z(self):
        raise NotImplementedError(str(type(self)) + " don't have implemented the function 'read_axis_z()'.")

    def read_axis_t(self, timestamp):
        raise NotImplementedError(str(type(self)) + " don't have implemented the function 'read_axis_t()'.")

    # Variables
    def read_variable_longitude(self):
        raise NotImplementedError(str(type(self)) + " don't have implemented the function 'read_variable_longitude()'.")

    def read_variable_latitude(self):
        raise NotImplementedError(str(type(self)) + " don't have implemented the function 'read_variable_latitude()'.")

    def read_variable_depth(self):
        raise NotImplementedError(str(type(self)) + " don't have implemented the function 'read_variable_depth()'.")

    def read_variable_time(self):
        raise NotImplementedError(str(type(self)) + " don't have implemented the function 'read_variable_time()'.")

    def read_variable_2D_sea_binary_mask(self):
        raise NotImplementedError(str(type(self)) + " don't have implemented the function 'read_variable_2D_sea_binary_mask()'.")

    def read_variable_2D_land_binary_mask(self):
        raise NotImplementedError(
            str(type(self)) + " don't have implemented the function 'read_variable_2D_land_binary_mask()'.")

    def read_variable_3D_sea_binary_mask(self):
        raise NotImplementedError(str(type(self)) + " don't have implemented the function 'read_variable_3D_sea_binary_mask()'.")

    def read_variable_3D_sea_binary_mask_at_time(self,index_t):
        raise NotImplementedError(str(type(self)) + " don't have implemented the function 'read_variable_3D_sea_binary_mask_at_time()'.")

    def read_variable_3D_land_binary_mask(self):
        raise NotImplementedError( str(type(self)) + " don't have implemented the function 'read_variable_3D_land_binary_mask()'.")

    def read_variable_3D_land_binary_mask_at_time(self,index_t):
        raise NotImplementedError(str(type(self)) + " don't have implemented the function 'read_variable_3D_land_binary_mas_at_time()'.")

    def read_variable_wet_binary_mask_at_time(self,index_t):
        raise NotImplementedError(str(type(self)) + " don't have implemented the function 'read_variable_wet_binary_mask_at_time()'.")

    def read_variable_mesh_size(self):
        raise NotImplementedError(str(type(self)) + " don't have implemented the function 'read_variable_mesh_size()'.")

    #################
    # HYDRO
    # Sea Surface
    #################
    def read_variable_sea_surface_height_above_mean_sea_level_at_time(self,index_t):
        raise NotImplementedError(str(type(self)) + " don't have implemented the function 'read_variable_sea_surface_height_above_mean_sea_level_at_time()'.")

    def read_variable_sea_surface_height_above_geoid_at_time(self,index_t):
        raise NotImplementedError(str(type(self)) + " don't have implemented the function 'read_variable_sea_surface_height_above_geoid_at_time()'.")

    def read_variable_sea_surface_temperature_at_time(self,index_t):
        raise NotImplementedError(str(type(self)) + " don't have implemented the function 'read_variable_sea_surface_temperature_at_time()'.")

    def read_variable_sea_surface_salinity_at_time(self,index_t):
        raise NotImplementedError(str(type(self)) + " don't have implemented the function 'read_variable_sea_surface_salinity_at_time()'.")

    def read_variable_sea_water_pressure_at_sea_water_surface_at_time(self,index_t):
        raise NotImplementedError(str(type(self)) + " don't have implemented the function 'read_variable_sea_water_pressure_at_sea_water_surface_at_time()'.")

    def read_variable_sea_surface_density_at_time(self,index_t):
        raise NotImplementedError(str(type(self)) + " don't have implemented the function 'read_variable_sea_surface_density_at_time()'.")

    def read_variable_sea_water_velocity_at_sea_water_surface_at_time(self,index_t):
        raise NotImplementedError(str(type(self)) + " don't have implemented the function 'read_variable_sea_water_velocity_at_sea_water_surface_at_time()'.")

    #################
    # HYDRO
    # Ground level
    #################
    def read_variable_sea_water_temperature_at_ground_level_at_time(self,index_t):
        raise NotImplementedError(str(type(self)) + " don't have implemented the function 'read_variable_sea_water_temperature_at_ground_level_at_time()'.")

    def read_variable_sea_water_salinity_at_ground_level_at_time(self,index_t):
        raise NotImplementedError(str(type(self)) + " don't have implemented the function 'read_variable_sea_water_salinity_at_ground_level_at_time()'.")

    def read_variable_sea_water_pressure_at_ground_level_at_time(self,index_t):
        raise NotImplementedError(str(type(self)) + " don't have implemented the function 'read_variable_sea_water_pressure_at_ground_level_at_time()'.")

    def read_variable_sea_water_density_at_ground_level_at_time(self,index_t):
        raise NotImplementedError(str(type(self)) + " don't have implemented the function 'read_variable_sea_water_density_at_ground_level_at_time()'.")

    def read_variable_sea_water_velocity_at_ground_level_at_time(self,index_t):
        raise NotImplementedError(str(type(self)) + " don't have implemented the function 'read_variable_sea_water_velocity_at_ground_level_at_time()'.")

    #################
    # HYDRO
    # 2D
    #################
    def read_variable_bathymetry(self):
        raise NotImplementedError(str(type(self)) + " don't have implemented the function 'read_variable_bathymetry()'.")

    def read_variable_barotropic_sea_water_velocity_at_time(self,index_t):
        raise NotImplementedError(str(type(self)) + " don't have implemented the function 'read_variable_barotropic_sea_water_velocity_at_time()'.")

    def read_variable_water_volume_transport_into_sea_water_from_rivers_at_time(self,index_t):
        raise NotImplementedError(str(type(self)) + " don't have implemented the function 'read_variable_water_volume_transport_into_sea_water_from_rivers_at_time()'.")

    #################
    # HYDRO
    # 3D
    #################
    def read_variable_sea_water_turbidity_at_time_and_depth(self,index_t,index_z):
        raise NotImplementedError(str(type(self)) + " don't have implemented the function 'read_variable_sea_water_turbidity_at_time_and_depth()'.")

    def read_variable_sea_water_electrical_conductivity_at_time_and_depth(self,index_t,index_z):
        raise NotImplementedError(str(type(self)) + " don't have implemented the function 'read_variable_sea_water_electrical_conductivity_at_time_and_depth()'.")

    def read_variable_sea_water_temperature_at_time_and_depth(self,index_t,index_z):
        raise NotImplementedError(str(type(self)) + " don't have implemented the function 'read_variable_sea_water_temperature_at_time_and_depth()'.")

    def read_variable_sea_water_salinity_at_time_and_depth(self,index_t,index_z):
        raise NotImplementedError(str(type(self)) + " don't have implemented the function 'read_variable_sea_water_salinity_at_time_and_depth()'.")

    def read_variable_sea_water_density_at_time_and_depth(self,index_t,index_z):
        raise NotImplementedError(str(type(self)) + " don't have implemented the function 'read_variable_sea_water_density_at_time_and_depth()'.")

    def read_variable_baroclinic_sea_water_velocity_at_time_and_depth(self,index_t,index_z):
        raise NotImplementedError(str(type(self)) + " don't have implemented the function 'read_variable_baroclinic_sea_water_velocity_at_time_and_depth()'.")

    #################
    # WAVES
    # Sea Surface
    #################
    def read_variable_sea_surface_wave_significant_height_at_time(self,index_t):
        raise NotImplementedError(str(type(self)) + " don't have implemented the function 'read_variable_sea_surface_wave_significant_height_at_time()'.")

    def read_variable_sea_surface_wave_breaking_height_at_time(self,index_t):
        raise NotImplementedError(str(type(self)) + " don't have implemented the function 'read_variable_sea_surface_wave_breaking_height_at_time()'.")

    def read_variable_sea_surface_wave_mean_period_at_time(self,index_t):
        raise NotImplementedError(str(type(self)) + " don't have implemented the function 'read_variable_sea_surface_wave_mean_period_at_time()'.")

    def read_variable_sea_surface_wave_peak_period_at_time(self,index_t):
        raise NotImplementedError(str(type(self)) + " don't have implemented the function 'read_variable_sea_surface_wave_peak_period_at_time()'.")

    def read_variable_sea_surface_wave_from_direction_at_time(self,index_t):
        raise NotImplementedError(str(type(self)) + " don't have implemented the function 'read_variable_sea_surface_wave_from_direction_at_time()'.")

    def read_variable_sea_surface_wave_to_direction_at_time(self,index_t):
        raise NotImplementedError(str(type(self)) + " don't have implemented the function 'read_variable_sea_surface_wave_to_direction_at_time()'.")

    def read_variable_sea_surface_wave_stokes_drift_velocity_at_time(self,index_t):
        raise NotImplementedError(str(type(self)) + " don't have implemented the function 'read_variable_sea_surface_wave_stokes_drift_velocity_at_time()'.")

    def read_variable_radiation_pressure_bernouilli_head_at_time(self,index_t):
        raise NotImplementedError(str(type(self)) + " don't have implemented the function 'read_variable_radiation_pressure_bernouilli_head_at_time()'.")

    def read_variable_sea_surface_wave_energy_flux_to_ocean_at_time(self,index_t):
        raise NotImplementedError(str(type(self)) + " don't have implemented the function 'read_variable_sea_surface_wave_energy_flux_to_ocean_at_time()'.")

    def read_variable_sea_surface_wave_energy_dissipation_at_ground_level_at_time(self,index_t):
        raise NotImplementedError(str(type(self)) + " don't have implemented the function 'read_variable_sea_surface_wave_energy_dissipation_at_ground_level_at_time()'.")

    #################
    # WAVES
    # Momentum flux
    #################
    def read_variable_atmosphere_momentum_flux_to_waves_at_time(self,index_t):
        raise NotImplementedError(str(type(self)) + " don't have implemented the function 'read_variable_atmosphere_momentum_flux_to_waves_at_time()'.")

    def read_variable_waves_momentum_flux_to_ocean_at_time(self,index_t):
        raise NotImplementedError(str(type(self)) + " don't have implemented the function 'read_variable_waves_momentum_flux_to_ocean_at_time()'.")

    #################
    # METEO
    # 2D
    #################
    def read_variable_topography(self):
        raise NotImplementedError(str(type(self)) + " don't have implemented the function 'read_variable_topography()'.")

    def read_variable_rainfall_amount_at_time(self,index_t):
        raise NotImplementedError(str(type(self)) + " don't have implemented the function 'read_variable_rainfall_amount_at_time()'.")

    #################
    # METEO
    # Surface air
    #################
    def read_variable_surface_air_pressure_at_time(self,index_t):
        raise NotImplementedError(str(type(self)) + " don't have implemented the function 'read_variable_surface_air_pressure_at_time()'.")

    def read_variable_sea_surface_air_pressure_at_time(self, index_t):
        raise NotImplementedError(str(type(self)) + " don't have implemented the function 'read_variable_sea_surface_air_pressure_at_time()'.")

    def read_variable_wind_stress_at_time(self,index_t):
        raise NotImplementedError(str(type(self)) + " don't have implemented the function 'read_variable_wind_stress_at_time()'.")

    def read_variable_surface_downward_sensible_heat_flux_at_time(self, index_t):
         raise NotImplementedError(str(type(self)) + " don't have implemented the function 'read_variable_surface_downward_sensible_heat_flux_at_time()'.")

    def read_variable_surface_downward_latent_heat_flux_at_time(self, index_t):
         raise NotImplementedError(str(type(self)) + " don't have implemented the function 'read_variable_surface_downward_latent_heat_flux_at_time()'.")

    def read_variable_surface_air_temperature_at_time(self, index_t):
         raise NotImplementedError(str(type(self)) + " don't have implemented the function 'read_variable_surface_air_temperature_at_time()'.")

    def read_variable_dew_point_temperature_at_time(self, index_t):
         raise NotImplementedError(str(type(self)) + " don't have implemented the function 'read_variable_dew_point_temperature_at_time()'.")

    def read_variable_surface_downwards_solar_radiation_at_time(self, index_t):
         raise NotImplementedError(str(type(self)) + " don't have implemented the function 'read_variable_surface_downwards_solar_radiation_at_time()'.")

    def read_variable_surface_downwards_thermal_radiation_at_time(self, index_t):
         raise NotImplementedError(str(type(self)) + " don't have implemented the function 'read_variable_surface_downwards_thermal_radiation_at_time()'.")

    def read_variable_surface_solar_radiation_at_time(self, index_t):
         raise NotImplementedError(str(type(self)) + " don't have implemented the function 'read_variable_surface_solar_radiation_at_time()'.")

    def read_variable_surface_thermal_radiation_at_time(self, index_t):
         raise NotImplementedError(str(type(self)) + " don't have implemented the function 'read_variable_surface_thermal_radiation_at_time()'.")

    #################
    # METEO
    # At 10 m
    #################
    def read_variable_wind_10m_at_time(self,index_t):
        raise NotImplementedError(str(type(self)) + " don't have implemented the function 'read_variable_eastward_wind_10m_at_time()'.")


