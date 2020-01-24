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
class CoverageWriter(object):
    
    def __init__(self, cov,myFile):
        self.coverage = cov;
        self.filename = myFile;

    def close(self):
        raise NotImplementedError(str(type(self)) + " don't have implemented the function 'close()'.")

    def write_variable_longitude(self):
        raise NotImplementedError(str(type(self)) + " don't have implemented the function 'write_variable_longitude()'.")

    def write_variable_latitude(self):
        raise NotImplementedError(str(type(self)) + " don't have implemented the function 'write_variable_latitude()'.")

    def write_variable_depth(self):
        raise NotImplementedError(str(type(self)) + " don't have implemented the function 'write_variable_depth()'.")

    def write_variable_time(self):
        raise NotImplementedError(str(type(self)) + " don't have implemented the function 'write_variable_time()'.")

    def write_variable_2D_sea_binary_mask(self, xmin, xmax, ymin, ymax):
        raise NotImplementedError(
            str(type(self)) + " don't have implemented the function 'write_variable_2D_sea_binary_mask()'.")

    def write_variable_2D_land_binary_mask(self, xmin, xmax, ymin, ymax):
        raise NotImplementedError(
            str(type(self)) + " don't have implemented the function 'write_variable_2D_land_binary_mask()'.")

    def write_variable_3D_sea_binary_mask(self, xmin, xmax, ymin, ymax):
        raise NotImplementedError(
            str(type(self)) + " don't have implemented the function 'write_variable_3D_sea_binary_mask()'.")

    def write_variable_3D_land_binary_mask(self, xmin, xmax, ymin, ymax):
        raise NotImplementedError(
            str(type(self)) + " don't have implemented the function 'write_variable_3D_land_binary_mask()'.")

    def write_variable_3D_land_binary_mask(self, xmin, xmax, ymin, ymax):
        raise NotImplementedError(
            str(type(self)) + " don't have implemented the function 'write_variable_3D_land_binary_mas()'.")

    def write_variable_wet_binary_mask(self, xmin, xmax, ymin, ymax):
        raise NotImplementedError(
            str(type(self)) + " don't have implemented the function 'write_variable_wet_binary_mask()'.")

    def write_variable_mesh_size(self, xmin, xmax, ymin, ymax):
        raise NotImplementedError(str(type(self)) + " don't have implemented the function 'write_variable_mesh_size()'.")

    #################
    # HYDRO
    # 2D
    #################
    def write_variable_bathymetry(self, xmin, xmax, ymin, ymax):
        raise NotImplementedError(
            str(type(self)) + " don't have implemented the function 'write_variable_bathymetry()'.")

    def write_variable_barotropic_sea_water_velocity(self):
        raise NotImplementedError(str(type(
            self)) + " don't have implemented the function 'write_variable_barotropic_sea_water_velocity()'.")

    def write_variable_water_volume_transport_into_sea_water_from_rivers(self):
        raise NotImplementedError(str(type(
            self)) + " don't have implemented the function 'write_variable_water_volume_transport_into_sea_water_from_rivers()'.")

    #################
    # HYDRO
    # Sea Surface
    #################
    def write_variable_sea_surface_height_above_mean_sea_level(self):
        raise NotImplementedError(str(type(
            self)) + " don't have implemented the function 'write_variable_sea_surface_height_above_mean_sea_level()'.")

    def write_variable_sea_surface_height_above_geoid(self):
        raise NotImplementedError(str(type(
            self)) + " don't have implemented the function 'write_variable_sea_surface_height_above_geoid()'.")

    def write_variable_sea_surface_temperature(self):
        raise NotImplementedError(
            str(type(self)) + " don't have implemented the function 'write_variable_sea_surface_temperature()'.")

    def write_variable_sea_surface_salinity(self):
        raise NotImplementedError(
            str(type(self)) + " don't have implemented the function 'write_variable_sea_surface_salinity()'.")

    def write_variable_sea_water_pressure_at_sea_water_surface(self):
        raise NotImplementedError(str(type(
            self)) + " don't have implemented the function 'write_variable_sea_water_pressure_at_sea_water_surface()'.")

    def write_variable_sea_surface_density(self):
        raise NotImplementedError(
            str(type(self)) + " don't have implemented the function 'write_variable_sea_surface_density()'.")

    def write_variable_sea_water_velocity_at_sea_water_surface(self):
        raise NotImplementedError(str(type(
            self)) + " don't have implemented the function 'write_variable_sea_water_velocity_at_sea_water_surface()'.")

    #################
    # HYDRO
    # Ground level
    #################
    def write_variable_sea_water_temperature_at_ground_level(self):
        raise NotImplementedError(str(type(
            self)) + " don't have implemented the function 'write_variable_sea_water_temperature_at_ground_level()'.")

    def write_variable_sea_water_salinity_at_ground_level(self):
        raise NotImplementedError(str(type(
            self)) + " don't have implemented the function 'write_variable_sea_water_salinity_at_ground_level()'.")

    def write_variable_sea_water_pressure_at_ground_level(self):
        raise NotImplementedError(str(type(
            self)) + " don't have implemented the function 'write_variable_sea_water_pressure_at_ground_level()'.")

    def write_variable_sea_water_density_at_ground_level(self):
        raise NotImplementedError(str(type(
            self)) + " don't have implemented the function 'write_variable_sea_water_density_at_ground_level()'.")

    def write_variable_sea_water_velocity_at_ground_level(self):
        raise NotImplementedError(str(type(
            self)) + " don't have implemented the function 'write_variable_sea_water_velocity_at_ground_level()'.")

    #################
    # HYDRO
    # 3D
    #################
    def write_variable_sea_water_turbidity_and_depth(self, index_z):
        raise NotImplementedError(str(type(
            self)) + " don't have implemented the function 'write_variable_sea_water_turbidity_and_depth()'.")

    def write_variable_sea_water_electrical_conductivity_and_depth(self, index_z):
        raise NotImplementedError(str(type(
            self)) + " don't have implemented the function 'write_variable_sea_water_electrical_conductivity_and_depth()'.")

    def write_variable_sea_water_temperature_and_depth(self, index_z):
        raise NotImplementedError(str(type(
            self)) + " don't have implemented the function 'write_variable_sea_water_temperature_and_depth()'.")

    def write_variable_sea_water_salinity_and_depth(self, index_z):
        raise NotImplementedError(str(type(
            self)) + " don't have implemented the function 'write_variable_sea_water_salinity_and_depth()'.")

    def write_variable_sea_water_density_and_depth(self, index_z):
        raise NotImplementedError(str(
            type(self)) + " don't have implemented the function 'write_variable_sea_water_density_and_depth()'.")

    def write_variable_baroclinic_sea_water_velocity_and_depth(self, index_z):
        raise NotImplementedError(str(type(
            self)) + " don't have implemented the function 'write_variable_baroclinic_sea_water_velocity_and_depth()'.")

    #################
    # WAVES
    # Sea Surface
    #################
    def write_variable_sea_surface_wave_significant_height(self):
        raise NotImplementedError(str(type(
            self)) + " don't have implemented the function 'write_variable_sea_surface_wave_significant_height()'.")

    def write_variable_sea_surface_wave_breaking_height(self):
        raise NotImplementedError(str(type(
            self)) + " don't have implemented the function 'write_variable_sea_surface_wave_breaking_height()'.")

    def write_variable_sea_surface_wave_mean_period(self):
        raise NotImplementedError(str(type(
            self)) + " don't have implemented the function 'write_variable_sea_surface_wave_mean_period()'.")

    def write_variable_sea_surface_wave_peak_period(self):
        raise NotImplementedError(str(type(
            self)) + " don't have implemented the function 'write_variable_sea_surface_wave_peak_period()'.")

    def write_variable_sea_surface_wave_from_direction(self):
        raise NotImplementedError(str(type(
            self)) + " don't have implemented the function 'write_variable_sea_surface_wave_from_direction()'.")

    def write_variable_sea_surface_wave_to_direction(self):
        raise NotImplementedError(str(type(
            self)) + " don't have implemented the function 'write_variable_sea_surface_wave_to_direction()'.")

    def write_variable_sea_surface_wave_stokes_drift_velocity(self):
        raise NotImplementedError(str(type(
            self)) + " don't have implemented the function 'write_variable_sea_surface_wave_stokes_drift_velocity()'.")

    def write_variable_radiation_pressure_bernouilli_head(self):
        raise NotImplementedError(str(type(
            self)) + " don't have implemented the function 'write_variable_radiation_pressure_bernouilli_head()'.")

    def write_variable_sea_surface_wave_energy_flux_to_ocean(self):
        raise NotImplementedError(str(type(
            self)) + " don't have implemented the function 'write_variable_sea_surface_wave_energy_flux_to_ocean()'.")

    #################
    # WAVES
    # Ground level
    #################
    def write_variable_sea_surface_wave_energy_dissipation_at_ground_level(self):
        raise NotImplementedError(str(type(
            self)) + " don't have implemented the function 'write_variable_sea_surface_wave_energy_dissipation_at_ground_level()'.")

    #################
    # WAVES
    # Momentum flux
    #################
    def write_variable_atmosphere_momentum_flux_to_waves(self):
        raise NotImplementedError(str(type(
            self)) + " don't have implemented the function 'write_variable_atmosphere_momentum_flux_to_waves()'.")

    def write_variable_waves_momentum_flux_to_ocean(self):
        raise NotImplementedError(str(type(
            self)) + " don't have implemented the function 'write_variable_waves_momentum_flux_to_ocean()'.")

    #################
    # METEO
    # 2D
    #################
    def write_variable_topography(self):
        raise NotImplementedError(
            str(type(self)) + " don't have implemented the function 'write_variable_topography()'.")

    def write_variable_rainfall_amount(self):
        raise NotImplementedError(
            str(type(self)) + " don't have implemented the function 'write_variable_rainfall_amount()'.")

    #################
    # METEO
    # Surface air
    #################
    def write_variable_surface_air_pressure(self):
        raise NotImplementedError(
            str(type(self)) + " don't have implemented the function 'write_variable_surface_air_pressure()'.")

    def write_variable_sea_surface_air_pressure(self):
        raise NotImplementedError(str(
            type(self)) + " don't have implemented the function 'write_variable_sea_surface_air_pressure()'.")

    def write_variable_wind_stress(self):
        raise NotImplementedError(
            str(type(self)) + " don't have implemented the function 'write_variable_wind_stress()'.")

    def write_variable_surface_downward_sensible_heat_flux(self):
        raise NotImplementedError(str(type(
            self)) + " don't have implemented the function 'write_variable_surface_downward_sensible_heat_flux()'.")

    def write_variable_surface_downward_latent_heat_flux(self):
        raise NotImplementedError(str(type(
            self)) + " don't have implemented the function 'write_variable_surface_downward_latent_heat_flux()'.")

    def write_variable_surface_air_temperature(self):
        raise NotImplementedError(
            str(type(self)) + " don't have implemented the function 'write_variable_surface_air_temperature()'.")

    def write_variable_dew_point_temperature(self):
        raise NotImplementedError(
            str(type(self)) + " don't have implemented the function 'write_variable_dew_point_temperature()'.")

    def write_variable_surface_downwards_solar_radiation(self):
        raise NotImplementedError(str(type(
            self)) + " don't have implemented the function 'write_variable_surface_downwards_solar_radiation()'.")

    def write_variable_surface_downwards_thermal_radiation(self):
        raise NotImplementedError(str(type(
            self)) + " don't have implemented the function 'write_variable_surface_downwards_thermal_radiation()'.")

    def write_variable_surface_solar_radiation(self):
        raise NotImplementedError(
            str(type(self)) + " don't have implemented the function 'write_variable_surface_solar_radiation()'.")

    def write_variable_surface_thermal_radiation(self):
        raise NotImplementedError(str(
            type(self)) + " don't have implemented the function 'write_variable_surface_thermal_radiation()'.")

    #################
    # METEO
    # At 10 m
    #################
    def write_variable_wind_10m(self):
        raise NotImplementedError(
            str(type(self)) + " don't have implemented the function 'write_variable_wind_10m()'.")

    def write_variable_wind_speed_10m(self):
        raise NotImplementedError(
            str(type(self)) + " don't have implemented the function 'write_variable_wind_speed_10m()'.")

    def write_variable_wind_to_direction_10m(self):
        raise NotImplementedError(
            str(type(self)) + " don't have implemented the function 'write_variable_wind_to_direction_10m()'.")

    def write_variable_wind_from_direction_10m(self):
        raise NotImplementedError(
            str(type(self)) + " don't have implemented the function 'write_variable_wind_from_direction_10m()'.")



