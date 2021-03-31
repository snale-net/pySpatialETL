from unittest import TestCase

from spatialetl.coverage import TimeCoverage
from spatialetl.coverage.io.netcdf.symphonie import SYMPHONIEReader


class TestTimeCoverage(TestCase):

    def setUp(self):
        reader = SYMPHONIEReader("resources/symphonie_grid.nc", "resources/symphonie_graphique.nc")
        self.coverage = TimeCoverage(reader)

    def test_init(self):
        self.assertIsNotNone(self.coverage)

    def test_create_mpi_map(self):
        self.fail()

    def test_update_mpi_map(self):
        self.fail()

    def test_find_time_index(self):
        self.fail()

    def test_read_axis_t(self):
        self.fail()

    def test_get_t_size(self):
        self.fail()

    def test_read_variable_2d_sea_binary_mask_at_time(self):
        self.fail()

    def test_read_variable_2d_wet_binary_mask_at_time(self):
        self.fail()

    def test_read_variable_2d_land_binary_mask_at_time(self):
        self.fail()

    def test_read_variable_sea_surface_height_above_mean_sea_level_at_time(self):
        self.fail()

    def test_read_variable_sea_surface_height_above_geoid_at_time(self):
        self.fail()

    def test_read_variable_sea_water_column_thickness_at_time(self):
        self.fail()

    def test_read_variable_sea_surface_temperature_at_time(self):
        self.fail()

    def test_read_variable_sea_surface_salinity_at_time(self):
        self.fail()

    def test_read_variable_sea_surface_pressure_at_time(self):
        self.fail()

    def test_read_variable_sea_surface_density_at_time(self):
        self.fail()

    def test_read_variable_sea_water_turbidity_at_time(self):
        self.fail()

    def test_read_variable_sea_water_velocity_at_sea_water_surface_at_time(self):
        self.fail()

    def test_read_variable_sea_water_temperature_at_ground_level_at_time(self):
        self.fail()

    def test_read_variable_sea_water_salinity_at_ground_level_at_time(self):
        self.fail()

    def test_read_variable_sea_water_velocity_at_ground_level_at_time(self):
        self.fail()

    def test_read_variable_barotropic_sea_water_velocity_at_time(self):
        self.fail()

    def test_read_variable_sea_surface_wave_significant_height_at_time(self):
        self.fail()

    def test_read_variable_sea_surface_wave_breaking_height_at_time(self):
        self.fail()

    def test_read_variable_sea_surface_wave_mean_period_at_time(self):
        self.fail()

    def test_read_variable_sea_surface_wave_peak_period_at_time(self):
        self.fail()

    def test_read_variable_sea_surface_wave_from_direction_at_time(self):
        self.fail()

    def test_read_variable_sea_surface_wave_to_direction_at_time(self):
        self.fail()

    def test_read_variable_sea_surface_wave_stokes_drift_velocity_at_time(self):
        self.fail()

    def test_read_variable_radiation_pressure_bernouilli_head_at_time(self):
        self.fail()

    def test_read_variable_sea_surface_wave_energy_flux_to_ocean_at_time(self):
        self.fail()

    def test_read_variable_sea_surface_wave_energy_dissipation_at_ground_level_at_time(self):
        self.fail()

    def test_read_variable_atmosphere_momentum_flux_to_waves_at_time(self):
        self.fail()

    def test_read_variable_waves_momentum_flux_to_ocean_at_time(self):
        self.fail()

    def test_read_variable_rainfall_amount_at_time(self):
        self.fail()

    def test_read_variable_surface_air_pressure_at_time(self):
        self.fail()

    def test_read_variable_sea_surface_air_pressure_at_time(self):
        self.fail()

    def test_read_variable_wind_stress_at_time(self):
        self.fail()

    def test_read_variable_surface_downward_sensible_heat_flux_at_time(self):
        self.fail()

    def test_read_variable_surface_downward_latent_heat_flux_at_time(self):
        self.fail()

    def test_read_variable_surface_air_temperature_at_time(self):
        self.fail()

    def test_read_variable_dew_point_temperature_at_time(self):
        self.fail()

    def test_read_variable_surface_downwards_solar_radiation_at_time(self):
        self.fail()

    def test_read_variable_surface_downwards_thermal_radiation_at_time(self):
        self.fail()

    def test_read_variable_surface_solar_radiation_at_time(self):
        self.fail()

    def test_read_variable_surface_thermal_radiation_at_time(self):
        self.fail()

    def test_read_variable_wind_10m_at_time(self):
        self.fail()

    def test_read_variable_wind_speed_10m_at_time(self):
        self.fail()

    def test_read_variable_wind_from_direction_10m_at_time(self):
        self.fail()

    def test_read_variable_wind_to_direction_10m_at_time(self):
        self.fail()
