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
# -
# Author : Fabien Rétif - fabien.retif@zoho.com
#
import sys

sys.path = ['/work/sciences/pySpatialETL'] + sys.path

from unittest import TestCase

import cftime
import numpy as np

from spatialetl.coverage.io.netcdf.symphonie.v293 import SYMPHONIEReader

# Lien vers le dossier de la lib



class TestSYMPHONIEReader(TestCase):

    def setUp(self):
        self.reader = SYMPHONIEReader("resources/grid.nc", "resources/2014*")

    def close(self):
        self.reader.close()

    def test_is_regular_grid(self):
        expected = False
        value = self.reader.is_regular_grid()
        self.assertEqual(expected,value)

    def test_get_x_size(self):
        expected = 12
        value = self.reader.get_x_size()
        self.assertEqual(expected,value)

    def test_get_y_size(self):
        expected = 12
        value = self.reader.get_y_size()
        self.assertEqual(expected,value)

    def test_get_z_size(self):
        expected = 10
        value = self.reader.get_z_size()
        self.assertEqual(expected,value)

    def test_get_t_size(self):
        expected = 3
        value = self.reader.get_t_size()
        self.assertEqual(expected,value)

    def test_read_axis_x(self):
        expected_shape = (5,5)
        expected_value = np.array([
        [-0.01798658, 0., 0.01798658, 0.03597315, 0.05395973],
        [-0.01798658, 0., 0.01798658, 0.03597315, 0.05395973],
        [-0.01798658, 0., 0.01798658, 0.03597315, 0.05395973],
        [-0.01798658, 0., 0.01798658, 0.03597315, 0.05395973],
        [-0.01798658, 0., 0.01798658, 0.03597315, 0.05395973]])

        candidate_value = self.reader.read_axis_x(0,5,0,5)
        self.assertEqual(expected_shape,np.shape(candidate_value),"shape assert")
        np.testing.assert_almost_equal(expected_value,candidate_value,5)

    def test_read_axis_y(self):
        expected_shape = (5, 5)
        expected_value = np.array([
        [-0.01798658, -0.01798658, -0.01798658, -0.01798658, -0.01798658],
        [ 0.,         0.,          0.,          0.,         0.        ],
        [ 0.01798658, 0.01798658, 0.01798658,  0.01798658,  0.01798658],
        [ 0.03597315, 0.03597315, 0.03597315, 0.03597315,  0.03597315],
        [ 0.05395972, 0.05395972, 0.05395972, 0.05395972,  0.05395972]])

        candidate_value = self.reader.read_axis_y(0, 5, 0, 5)
        self.assertEqual(expected_shape, np.shape(candidate_value), "shape assert")
        np.testing.assert_almost_equal(expected_value, candidate_value, 5)

    def test_read_axis_z(self):
        expected_shape = (10, 12, 12)
        expected_value = np.array([
            [-0.01798658, -0.01798658, -0.01798658, -0.01798658, -0.01798658],
            [0., 0., 0., 0., 0.],
            [0.01798658, 0.01798658, 0.01798658, 0.01798658, 0.01798658],
            [0.03597315, 0.03597315, 0.03597315, 0.03597315, 0.03597315],
            [0.05395972, 0.05395972, 0.05395972, 0.05395972, 0.05395972]])

        candidate_value = self.reader.read_axis_z()
        self.assertEqual(expected_shape, np.shape(candidate_value), "shape assert")
        #np.testing.assert_almost_equal(expected_value, candidate_value, 5)

    def test_read_axis_t(self):
        expected_shape = (2,)
        expected_value = np.array([1388579767.0, 1388581212.0])

        candidate_value = self.reader.read_axis_t(1,3,timestamp=1)
        self.assertEqual(expected_shape, np.shape(candidate_value), "shape assert")
        np.testing.assert_almost_equal(expected_value, candidate_value, 5)

    # Variables
    def test_read_variable_longitude(self):
        expected_shape = (12, 12)
        expected_value = np.array([
        [-0.01798658,  0.,          0.01798658,  0.03597315,  0.05395973,  0.0719463,
        0.08993288,  0.10791946,  0.12590603,  0.14389261,  0.16187918,  0.17986576],
        [-0.01798658,  0.,          0.01798658,  0.03597315,  0.05395973,  0.0719463,
        0.08993288,  0.10791946,  0.12590603,  0.14389261,  0.16187918,  0.17986576],
        [-0.01798658,  0.,          0.01798658,  0.03597315,  0.05395973,  0.0719463,
        0.08993288,  0.10791946,  0.12590603,  0.14389261,  0.16187918,  0.17986576],
        [-0.01798658,  0.,          0.01798658,  0.03597315,  0.05395973,  0.0719463,
        0.08993288,  0.10791946,  0.12590603,  0.14389261,  0.16187918,  0.17986576],
        [-0.01798658,  0.,          0.01798658,  0.03597315,  0.05395973,  0.0719463,
        0.08993288,  0.10791946,  0.12590603,  0.14389261,  0.16187918,  0.17986576],
        [-0.01798658,  0.,          0.01798658,  0.03597315,  0.05395973,  0.0719463,
        0.08993288,  0.10791946,  0.12590603,  0.14389261,  0.16187918,  0.17986576],
        [-0.01798658,  0.,          0.01798658,  0.03597315,  0.05395973,  0.0719463,
        0.08993288,  0.10791946,  0.12590603,  0.14389261,  0.16187918,  0.17986576],
        [-0.01798658,  0.,          0.01798658,  0.03597315,  0.05395973,  0.0719463,
        0.08993288,  0.10791946,  0.12590603,  0.14389261,  0.16187918,  0.17986576],
        [-0.01798658,  0.,          0.01798658,  0.03597315,  0.05395973,  0.0719463,
        0.08993288,  0.10791946,  0.12590603,  0.14389261,  0.16187918,  0.17986576],
        [-0.01798658,  0.,          0.01798658,  0.03597315,  0.05395973,  0.0719463,
        0.08993288,  0.10791946,  0.12590603,  0.14389261,  0.16187918,  0.17986576],
        [-0.01798658,  0.,          0.01798658,  0.03597315,  0.05395973,  0.0719463,
        0.08993288,  0.10791946,  0.12590603,  0.14389261, 0.16187918,  0.17986576],
        [-0.01798658,  0.,          0.01798658,  0.03597315,  0.05395973,  0.0719463,
        0.08993288,  0.10791946, 0.12590603,  0.14389261,  0.16187918,  0.17986576]])

        candidate_value = self.reader.read_variable_longitude(0,12,0,12)
        self.assertEqual(expected_shape, np.shape(candidate_value), "shape assert")
        np.testing.assert_almost_equal(expected_value, candidate_value, 5)

    def test_read_variable_latitude(self):
        expected_shape = (2, 2)
        expected_value = np.array([[0.07194629, 0.07194629],[0.08993284, 0.08993284]])
        candidate_value = self.reader.read_variable_latitude(5,7,5,7)
        self.assertEqual(expected_shape, np.shape(candidate_value), "shape assert")
        np.testing.assert_almost_equal(expected_value, candidate_value, 5)

    def test_read_variable_time(self):
        expected_shape = (2,)
        expected_value = np.array([cftime.datetime(2014,1,1,12,36,7), cftime.datetime(2014,1,1,13,0,12)])

        candidate_value = self.reader.read_variable_time(1,3,timestamp=0)
        self.assertEqual(expected_shape, np.shape(candidate_value), "shape assert")

        for i in range(0, len(expected_value)):
            self.assertEqual(expected_value[i],candidate_value[i], "datetime assert")

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

