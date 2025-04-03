spatialetl.coverage.io.CoverageWriter
=====================================

.. py:module:: spatialetl.coverage.io.CoverageWriter


Classes
-------

.. autoapisummary::

   spatialetl.coverage.io.CoverageWriter.CoverageWriter


Module Contents
---------------

.. py:class:: CoverageWriter(cov, myFile)

   Bases: :py:obj:`object`


   .. py:attribute:: coverage


   .. py:attribute:: filename


   .. py:method:: close()
      :abstractmethod:



   .. py:method:: write_variable_longitude()
      :abstractmethod:



   .. py:method:: write_variable_latitude()
      :abstractmethod:



   .. py:method:: write_variable_depth()
      :abstractmethod:



   .. py:method:: write_variable_time()
      :abstractmethod:



   .. py:method:: write_variable_2D_sea_binary_mask()
      :abstractmethod:



   .. py:method:: write_variable_2D_land_binary_mask()
      :abstractmethod:



   .. py:method:: write_variable_3D_sea_binary_mask()
      :abstractmethod:



   .. py:method:: write_variable_3D_land_binary_mask()
      :abstractmethod:



   .. py:method:: write_variable_wet_binary_mask()
      :abstractmethod:



   .. py:method:: write_variable_mesh_size()
      :abstractmethod:



   .. py:method:: write_variable_bathymetry()
      :abstractmethod:



   .. py:method:: write_variable_barotropic_sea_water_velocity()
      :abstractmethod:



   .. py:method:: write_variable_water_volume_transport_into_sea_water_from_rivers()
      :abstractmethod:



   .. py:method:: write_variable_sea_surface_height_above_mean_sea_level()
      :abstractmethod:



   .. py:method:: write_variable_sea_surface_height_above_geoid()
      :abstractmethod:



   .. py:method:: write_variable_sea_surface_temperature()
      :abstractmethod:



   .. py:method:: write_variable_sea_surface_salinity()
      :abstractmethod:



   .. py:method:: write_variable_sea_water_pressure_at_sea_water_surface()
      :abstractmethod:



   .. py:method:: write_variable_sea_surface_density()
      :abstractmethod:



   .. py:method:: write_variable_sea_water_velocity_at_sea_water_surface()
      :abstractmethod:



   .. py:method:: write_variable_sea_water_temperature_at_ground_level()
      :abstractmethod:



   .. py:method:: write_variable_sea_water_salinity_at_ground_level()
      :abstractmethod:



   .. py:method:: write_variable_sea_water_pressure_at_ground_level()
      :abstractmethod:



   .. py:method:: write_variable_sea_water_density_at_ground_level()
      :abstractmethod:



   .. py:method:: write_variable_sea_water_velocity_at_ground_level()
      :abstractmethod:



   .. py:method:: write_variable_sea_water_turbidity_and_depth(index_z)
      :abstractmethod:



   .. py:method:: write_variable_sea_water_electrical_conductivity_and_depth(index_z)
      :abstractmethod:



   .. py:method:: write_variable_sea_water_temperature_and_depth(index_z)
      :abstractmethod:



   .. py:method:: write_variable_sea_water_salinity_and_depth(index_z)
      :abstractmethod:



   .. py:method:: write_variable_sea_water_density_and_depth(index_z)
      :abstractmethod:



   .. py:method:: write_variable_baroclinic_sea_water_velocity_and_depth(index_z)
      :abstractmethod:



   .. py:method:: write_variable_sea_surface_wave_significant_height()
      :abstractmethod:



   .. py:method:: write_variable_sea_surface_wave_breaking_height()
      :abstractmethod:



   .. py:method:: write_variable_sea_surface_wave_mean_period()
      :abstractmethod:



   .. py:method:: write_variable_sea_surface_wave_peak_period()
      :abstractmethod:



   .. py:method:: write_variable_sea_surface_wave_from_direction()
      :abstractmethod:



   .. py:method:: write_variable_sea_surface_wave_to_direction()
      :abstractmethod:



   .. py:method:: write_variable_sea_surface_wave_stokes_drift_velocity()
      :abstractmethod:



   .. py:method:: write_variable_radiation_pressure_bernouilli_head()
      :abstractmethod:



   .. py:method:: write_variable_sea_surface_wave_energy_flux_to_ocean()
      :abstractmethod:



   .. py:method:: write_variable_sea_surface_wave_energy_dissipation_at_ground_level()
      :abstractmethod:



   .. py:method:: write_variable_atmosphere_momentum_flux_to_waves()
      :abstractmethod:



   .. py:method:: write_variable_waves_momentum_flux_to_ocean()
      :abstractmethod:



   .. py:method:: write_variable_topography()
      :abstractmethod:



   .. py:method:: write_variable_rainfall_amount()
      :abstractmethod:



   .. py:method:: write_variable_surface_air_pressure()
      :abstractmethod:



   .. py:method:: write_variable_sea_surface_air_pressure()
      :abstractmethod:



   .. py:method:: write_variable_wind_stress()
      :abstractmethod:



   .. py:method:: write_variable_surface_downward_sensible_heat_flux()
      :abstractmethod:



   .. py:method:: write_variable_surface_downward_latent_heat_flux()
      :abstractmethod:



   .. py:method:: write_variable_surface_air_temperature()
      :abstractmethod:



   .. py:method:: write_variable_dew_point_temperature()
      :abstractmethod:



   .. py:method:: write_variable_surface_downward_solar_radiation()
      :abstractmethod:



   .. py:method:: write_variable_surface_downward_thermal_radiation()
      :abstractmethod:



   .. py:method:: write_variable_surface_solar_radiation()
      :abstractmethod:



   .. py:method:: write_variable_surface_thermal_radiation()
      :abstractmethod:



   .. py:method:: write_variable_wind_10m()
      :abstractmethod:



   .. py:method:: write_variable_wind_speed_10m()
      :abstractmethod:



   .. py:method:: write_variable_wind_to_direction_10m()
      :abstractmethod:



   .. py:method:: write_variable_wind_from_direction_10m()
      :abstractmethod:



