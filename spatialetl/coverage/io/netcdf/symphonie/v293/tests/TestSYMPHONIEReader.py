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
from unittest import TestCase

import numpy as np

from spatialetl.coverage.io.netcdf.symphonie.v293 import SYMPHONIEReader


class SymphonieReaderTest(TestCase):

    def setUp(self):
        self.reader = SYMPHONIEReader("cropped_grid.nc", "cropped_example.nc")

    def close(self):
        self.reader.close()

    def test_is_regular_grid(self):
        expected = False
        value = self.reader.is_regular_grid()
        self.assertEqual(expected,value)

    def test_get_x_size(self):
        expected = 5
        value = self.reader.get_x_size()
        self.assertEqual(expected,value)

    def test_get_y_size(self):
        expected = 5
        value = self.reader.get_y_size()
        self.assertEqual(expected,value)

    def test_get_z_size(self):
        expected = 10
        value = self.reader.get_z_size()
        self.assertEqual(expected,value)

    def test_get_t_size(self):
        expected = 1
        value = self.reader.get_t_size()
        self.assertEqual(expected,value)

    def test_read_axis_x(self):
        expected_shape = (5,5)
        expected_value = np.array([
        [3.01484219,3.0152151,3.01558801,3.0159549,3.0163157],
        [3.0150673,3.01544609,3.01582489,3.01619746,3.01656373],
        [3.01529241,3.01567709,3.01606177,3.01644002,3.01681176],
        [3.01552402,3.01591467,3.01630532,3.01668933,3.01706662],
        [3.01576238,3.01615909,3.0165558,3.01694566,3.01732856]])

        value = self.reader.read_axis_x(0,5,0,5)
        self.assertEqual(expected_shape,np.shape(value),"shape assert")
        #self.assertEqual(expected_value.tolist(), value.tolist(),"value assert")

    def test_read_axis_y(self):
        expected_shape = (5, 5)
        value = self.reader.read_axis_y(0, 5, 0, 5)
        self.assertEqual(expected_shape, np.shape(value), "shape assert")

    def test_read_axis_z(self):
        expected_shape = (10, 5, 5)
        value = self.reader.read_axis_z()
        self.assertEqual(expected_shape, np.shape(value), "shape assert")

    def test_read_axis_t(self):
        expected_shape = (0,)
        value = self.reader.read_axis_t(0,1,timestamp=1)
        self.assertEqual(expected_shape, np.shape(value), "shape assert")

    # Variables
    def test_read_variable_longitude(self):
        value = self.reader.read_variable_longitude(0,5,0,5)

    def test_read_variable_latitude(self):
        value = self.reader.read_variable_latitude(0,5,0,5)

    def test_read_variable_time(self):
        self.reader.read_variable_time(0,1)

    def test_read_variable_2D_sea_binary_mask(self):
        expected_value = np.array([
            [12.0, 13.0],
            [17.0, 18.0],
        ])
        value = self.reader.read_variable_2D_wet_binary_mask_at_time(0, 2, 4, 2, 4);
        self.assertEqual(expected_value.tolist(), value.tolist())

    def test_read_variable_2D_land_binary_mask(self):
        value = self.reader.read_variable_2D_land_binary_mask_at_time(0, 2, 4, 2, 4);

    def test_read_variable_3D_sea_binary_mask(self):
        value = self.reader.read_variable_2D_sea_binary_mask_at_time(0, 2, 4, 2, 4);

    def test_read_variable_3D_sea_binary_mask_at_time(self):
        value = self.reader.read_variable_3D_sea_binary_mask_at_time(0, 2, 4, 2, 4);

    def test_read_variable_3D_land_binary_mask(self):
        value = self.reader.read_variable_3D_land_binary_mask(2, 4, 2, 4);

    def test_read_variable_3D_land_binary_mask_at_time(self):
        value = self.reader.read_variable_3D_land_binary_mask_at_time(0, 2, 4, 2, 4);

    def test_read_variable_2D_wet_binary_mask_at_time(self):
        expected = np.array([
            [12.0, 13.0],
            [17.0, 18.0],
        ])
        values = self.reader.read_variable_2D_wet_binary_mask_at_time(0, 2, 4, 2, 4);
        self.assertEqual(expected.tolist(), values.tolist())

    def test_read_variable_mesh_size(self):
        value = self.reader.read_variable_mesh_size(2, 4, 2, 4);

    #################
    # HYDRO
    # 2D
    #################
    def test_read_variable_bathymetry(self):
        value = self.reader.read_variable_bathymetry(2, 4, 2, 4);

    def test_read_variable_barotropic_sea_water_velocity_at_time(self):
        value = self.reader.read_variable_barotropic_sea_water_velocity_at_time(0,2, 4, 2, 4);

    def test_read_variable_water_volume_transport_into_sea_water_from_rivers_at_time(self):
        value = self.reader.read_variable_water_volume_transport_into_sea_water_from_rivers_at_time()

    #################
    # HYDRO
    # Sea Surface
    #################
    def test_read_variable_sea_surface_height_above_mean_sea_level_at_time(self):
        expected = np.array([
            [12.0, 13.0],
            [17.0, 18.0],
        ])
        value = self.reader.read_variable_sea_surface_height_above_mean_sea_level_at_time(0, 2, 4, 2, 4);
        self.assertEqual(expected.tolist(), value.tolist())

    def tes_read_variable_sea_water_column_thickness_at_time(self):
        expected = np.array([
            [12.0, 13.0],
            [17.0, 18.0],
        ])
        value = self.reader.read_variable_sea_water_column_thickness_at_time(0, 2, 4, 2, 4);
        self.assertEqual(expected.tolist(), value.tolist())

    def test_read_variable_sea_surface_height_above_geoid_at_time(self):
        value = self.reader.read_variable_sea_surface_height_above_geoid_at_time(0, 2, 4, 2, 4);

    def test_read_variable_sea_surface_temperature_at_time(self):
        expected = np.array([
            [12.0, 13.0],
            [17.0, 18.0],
        ])
        value = self.reader.read_variable_sea_surface_temperature_at_time(0, 2, 4, 2, 4);
        self.assertEqual(expected.tolist(), value.tolist())

    def test_read_variable_sea_surface_salinity_at_time(self):
        expected = np.array([
            [12.0, 13.0],
            [17.0, 18.0],
        ])
        value = self.reader.read_variable_sea_surface_salinity_at_time(0, 2, 4, 2, 4);
        self.assertEqual(expected.tolist(), value.tolist())

    def test_read_variable_sea_water_pressure_at_sea_water_surface_at_time(self):
        value = self.reader.read_variable_sea_water_pressure_at_sea_water_surface_at_time(0, 2, 4, 2, 4);

    def test_read_variable_sea_surface_density_at_time(self):
        value = self.reader.read_variable_sea_surface_density_at_time(0, 2, 4, 2, 4);

    def test_read_variable_sea_water_velocity_at_sea_water_surface_at_time(self):
       value = self.reader.read_variable_sea_water_velocity_at_sea_water_surface_at_time(0,0,5,0,5)

    #################
    # HYDRO
    # Ground level
    #################
    def test_read_variable_sea_water_temperature_at_ground_level_at_time(self):
        value = self.reader.read_variable_sea_water_temperature_at_ground_level_at_time(0, 0, 5, 0, 5)

    def test_read_variable_sea_water_salinity_at_ground_level_at_time(self):
        value = self.reader.read_variable_sea_water_salinity_at_ground_level_at_time(0, 0, 5, 0, 5)

    def test_read_variable_sea_water_pressure_at_ground_level_at_time(self):
        value = self.reader.read_variable_sea_water_pressure_at_ground_level_at_time(0, 0, 5, 0, 5)

    def test_read_variable_sea_water_density_at_ground_level_at_time(self):
        value = self.reader.read_variable_sea_water_density_at_ground_level_at_time(0, 0, 5, 0, 5)

    def test_read_variable_sea_water_velocity_at_ground_level_at_time(self):
        value = self.reader.read_variable_sea_water_velocity_at_ground_level_at_time(0, 0, 5, 0, 5)

    #################
    # HYDRO
    # 3D
    #################
    def test_read_variable_depth_at_depth(self):
        value = self.reader.read_variable_depth_at_depth(0, 5, 0, 5)

    def test_read_variable_sea_water_turbidity_at_time_and_depth(self):
        value = self.reader.read_variable_sea_water_turbidity_at_time_and_depth(0, 0, 0 ,5 ,0 ,5)

    def test_read_variable_sea_water_electrical_conductivity_at_time_and_depth(self):
        value = self.reader.read_variable_sea_water_electrical_conductivity_at_time_and_depth(0, 0, 0, 5, 0, 5)

    def test_read_variable_sea_water_temperature_at_time_and_depth(self):
        value = self.reader.read_variable_sea_water_temperature_at_time_and_depth(0, 0, 0, 5, 0, 5)

    def test_read_variable_sea_water_salinity_at_time_and_depth(self):
        value = self.reader.read_variable_sea_water_salinity_at_time_and_depth(0, 0, 0, 5, 0, 5)

    def test_read_variable_sea_water_density_at_time_and_depth(self):
        value = self.reader.read_variable_sea_water_density_at_time_and_depth(0, 0, 0, 5, 0, 5)

    def test_read_variable_baroclinic_sea_water_velocity_at_time_and_depth(self):
        value = self.reader.read_variable_baroclinic_sea_water_velocity_at_time_and_depth(0, 0, 0, 5, 0, 5)

    #################
    # WAVES
    # Sea Surface
    #################
    def test_read_variable_sea_surface_wave_significant_height_at_time(self):
        value = self.reader.read_variable_sea_surface_wave_significant_height_at_time(0, 0, 5, 0, 5)

    def test_read_variable_sea_surface_wave_breaking_height_at_time(self):
        value = self.reader.read_variable_sea_surface_wave_breaking_height_at_time(0, 0, 5, 0, 5)

    def test_read_variable_sea_surface_wave_mean_period_at_time(self):
        value = self.reader.read_variable_sea_surface_wave_mean_period_at_time(0, 0, 5, 0, 5)

    def test_read_variable_sea_surface_wave_peak_period_at_time(self):
        value = self.reader.read_variable_sea_surface_wave_peak_period_at_time(0, 0, 5, 0, 5)

    def test_read_variable_sea_surface_wave_from_direction_at_time(self):
        value = self.reader.read_variable_sea_surface_wave_from_direction_at_time(0, 0, 5, 0, 5)

    def test_read_variable_sea_surface_wave_to_direction_at_time(self):
        value = self.reader.read_variable_sea_surface_wave_to_direction_at_time(0, 0, 5, 0, 5)

    def test_read_variable_sea_surface_wave_stokes_drift_velocity_at_time(self):
        value = self.reader.read_variable_sea_surface_wave_stokes_drift_velocity_at_time(0, 0, 5, 0, 5)

    def test_read_variable_radiation_pressure_bernouilli_head_at_time(self):
        value = self.reader.read_variable_radiation_pressure_bernouilli_head_at_time(0, 0, 5, 0, 5)

    def test_read_variable_sea_surface_wave_energy_flux_to_ocean_at_time(self):
        value = self.reader.read_variable_sea_surface_wave_energy_flux_to_ocean_at_time(0, 0, 5, 0, 5)

    #################
    # WAVES
    # Ground level
    #################
    def test_read_variable_sea_surface_wave_energy_dissipation_at_ground_level_at_time(self):
        value = self.reader.read_variable_sea_surface_wave_energy_dissipation_at_ground_level_at_time(0, 0, 5, 0, 5)

    #################
    # WAVES
    # Momentum flux
    #################
    def test_read_variable_atmosphere_momentum_flux_to_waves_at_time(self):
        value = self.reader.read_variable_atmosphere_momentum_flux_to_waves_at_time(0, 0, 5, 0, 5)

    def test_read_variable_waves_momentum_flux_to_ocean_at_time(self):
        value = self.reader.read_variable_waves_momentum_flux_to_ocean_at_time(0, 0, 5, 0, 5)

    #################
    # METEO
    # 2D
    #################
    def test_read_variable_topography(self):
        value = self.reader.read_variable_topography(0, 5, 0, 5)

    def test_read_variable_rainfall_amount_at_time(self):
        value = self.reader.read_variable_rainfall_amount_at_time(0, 0, 5, 0, 5)

    #################
    # METEO
    # Surface air
    #################
    def test_read_variable_surface_air_pressure_at_time(self):
        value = self.reader.read_variable_surface_air_pressure_at_time(0, 0, 5, 0, 5)

    def test_read_variable_sea_surface_air_pressure_at_time(self):
        value = self.reader.read_variable_sea_surface_air_pressure_at_time(0, 0, 5, 0, 5)

    def test_read_variable_wind_stress_at_time(self):
        value = self.reader.read_variable_wind_stress_at_time(0, 0, 5, 0, 5)

    def test_read_variable_surface_downward_sensible_heat_flux_at_time(self):
        value = self.reader.read_variable_surface_downward_sensible_heat_flux_at_time(0, 0, 5, 0, 5)

    def test_read_variable_surface_downward_latent_heat_flux_at_time(self):
        value = self.reader.read_variable_surface_downward_latent_heat_flux_at_time(0, 0, 5, 0, 5)

    def test_read_variable_surface_air_temperature_at_time(self):
        value = self.reader.read_variable_surface_air_temperature_at_time(0, 0, 5, 0, 5)

    def test_read_variable_dew_point_temperature_at_time(self):
        value = self.reader.read_variable_dew_point_temperature_at_time(0, 0, 5, 0, 5)

    def test_read_variable_surface_downwards_solar_radiation_at_time(self):
        value = self.reader.read_variable_surface_downwards_solar_radiation_at_time(0, 0, 5, 0, 5)

    def test_read_variable_surface_downwards_thermal_radiation_at_time(self):
        value = self.reader.read_variable_surface_downwards_thermal_radiation_at_time(0, 0, 5, 0, 5)

    def test_read_variable_surface_solar_radiation_at_time(self):
        value = self.reader.read_variable_surface_solar_radiation_at_time(0, 0, 5, 0, 5)

    def test_read_variable_surface_thermal_radiation_at_time(self):
        value = self.reader.read_variable_surface_thermal_radiation_at_time(0, 0, 5, 0, 5)

    #################
    # METEO
    # At 10 m
    #################
    def test_read_variable_wind_10m_at_time(self):
        value = self.reader.read_variable_wind_10m_at_time(0, 0, 5, 0, 5)

