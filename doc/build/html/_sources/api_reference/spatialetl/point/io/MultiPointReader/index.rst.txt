spatialetl.point.io.MultiPointReader
====================================

.. py:module:: spatialetl.point.io.MultiPointReader


Classes
-------

.. autoapisummary::

   spatialetl.point.io.MultiPointReader.MultiPointReader


Module Contents
---------------

.. py:class:: MultiPointReader(myFile)

   Bases: :py:obj:`object`


   .. py:attribute:: filename


   .. py:method:: close()
      :abstractmethod:



   .. py:method:: get_x_size()
      :abstractmethod:



   .. py:method:: get_y_size()
      :abstractmethod:



   .. py:method:: get_z_size()
      :abstractmethod:



   .. py:method:: get_t_size()
      :abstractmethod:



   .. py:method:: read_axis_x()
      :abstractmethod:



   .. py:method:: read_axis_y()
      :abstractmethod:



   .. py:method:: read_axis_z()
      :abstractmethod:



   .. py:method:: read_axis_t(tmin, tmax, timestamp)
      :abstractmethod:



   .. py:method:: read_variable_point_names()


   .. py:method:: read_variable_longitude()
      :abstractmethod:



   .. py:method:: read_variable_latitude()
      :abstractmethod:



   .. py:method:: read_variable_depth()
      :abstractmethod:



   .. py:method:: read_variable_time()
      :abstractmethod:



   .. py:method:: read_variable_2D_sea_binary_mask()
      :abstractmethod:



   .. py:method:: read_variable_2D_land_binary_mask()
      :abstractmethod:



   .. py:method:: read_variable_3D_sea_binary_mask()
      :abstractmethod:



   .. py:method:: read_variable_3D_sea_binary_mask_at_time(index_t)
      :abstractmethod:



   .. py:method:: read_variable_3D_land_binary_mask()
      :abstractmethod:



   .. py:method:: read_variable_3D_land_binary_mask_at_time(index_t)
      :abstractmethod:



   .. py:method:: read_variable_wet_binary_mask_at_time(index_t)
      :abstractmethod:



   .. py:method:: read_variable_mesh_size()
      :abstractmethod:



   .. py:method:: read_variable_sea_surface_height_above_mean_sea_level_at_time(index_t)
      :abstractmethod:



   .. py:method:: read_variable_sea_surface_height_above_geoid_at_time(index_t)
      :abstractmethod:



   .. py:method:: read_variable_sea_surface_temperature_at_time(index_t)
      :abstractmethod:



   .. py:method:: read_variable_sea_surface_salinity_at_time(index_t)
      :abstractmethod:



   .. py:method:: read_variable_sea_water_pressure_at_sea_water_surface_at_time(index_t)
      :abstractmethod:



   .. py:method:: read_variable_sea_surface_density_at_time(index_t)
      :abstractmethod:



   .. py:method:: read_variable_sea_water_velocity_at_sea_water_surface_at_time(index_t)
      :abstractmethod:



   .. py:method:: read_variable_sea_water_temperature_at_ground_level_at_time(index_t)
      :abstractmethod:



   .. py:method:: read_variable_sea_water_salinity_at_ground_level_at_time(index_t)
      :abstractmethod:



   .. py:method:: read_variable_sea_water_pressure_at_ground_level_at_time(index_t)
      :abstractmethod:



   .. py:method:: read_variable_sea_water_density_at_ground_level_at_time(index_t)
      :abstractmethod:



   .. py:method:: read_variable_sea_water_velocity_at_ground_level_at_time(index_t)
      :abstractmethod:



   .. py:method:: read_variable_bathymetry()
      :abstractmethod:



   .. py:method:: read_variable_barotropic_sea_water_velocity_at_time(index_t)
      :abstractmethod:



   .. py:method:: read_variable_water_volume_transport_into_sea_water_from_rivers_at_time(index_t)
      :abstractmethod:



   .. py:method:: read_variable_sea_water_turbidity_at_time_and_depth(index_t, index_z)
      :abstractmethod:



   .. py:method:: read_variable_sea_water_electrical_conductivity_at_time_and_depth(index_t, index_z)
      :abstractmethod:



   .. py:method:: read_variable_sea_water_temperature_at_time_and_depth(index_t, index_z)
      :abstractmethod:



   .. py:method:: read_variable_sea_water_salinity_at_time_and_depth(index_t, index_z)
      :abstractmethod:



   .. py:method:: read_variable_sea_water_density_at_time_and_depth(index_t, index_z)
      :abstractmethod:



   .. py:method:: read_variable_baroclinic_sea_water_velocity_at_time_and_depth(index_t, index_z)
      :abstractmethod:



   .. py:method:: read_variable_sea_surface_wave_significant_height_at_time(index_t)
      :abstractmethod:



   .. py:method:: read_variable_sea_surface_wave_breaking_height_at_time(index_t)
      :abstractmethod:



   .. py:method:: read_variable_sea_surface_wave_mean_period_at_time(index_t)
      :abstractmethod:



   .. py:method:: read_variable_sea_surface_wave_peak_period_at_time(index_t)
      :abstractmethod:



   .. py:method:: read_variable_sea_surface_wave_from_direction_at_time(index_t)
      :abstractmethod:



   .. py:method:: read_variable_sea_surface_wave_to_direction_at_time(index_t)
      :abstractmethod:



   .. py:method:: read_variable_sea_surface_wave_stokes_drift_velocity_at_time(index_t)
      :abstractmethod:



   .. py:method:: read_variable_radiation_pressure_bernouilli_head_at_time(index_t)
      :abstractmethod:



   .. py:method:: read_variable_sea_surface_wave_energy_flux_to_ocean_at_time(index_t)
      :abstractmethod:



   .. py:method:: read_variable_sea_surface_wave_energy_dissipation_at_ground_level_at_time(index_t)
      :abstractmethod:



   .. py:method:: read_variable_atmosphere_momentum_flux_to_waves_at_time(index_t)
      :abstractmethod:



   .. py:method:: read_variable_waves_momentum_flux_to_ocean_at_time(index_t)
      :abstractmethod:



   .. py:method:: read_variable_topography()
      :abstractmethod:



   .. py:method:: read_variable_rainfall_amount_at_time(index_t)
      :abstractmethod:



   .. py:method:: read_variable_surface_air_pressure_at_time(index_t)
      :abstractmethod:



   .. py:method:: read_variable_sea_surface_air_pressure_at_time(index_t)
      :abstractmethod:



   .. py:method:: read_variable_wind_stress_at_time(index_t)
      :abstractmethod:



   .. py:method:: read_variable_surface_downward_sensible_heat_flux_at_time(index_t)
      :abstractmethod:



   .. py:method:: read_variable_surface_downward_latent_heat_flux_at_time(index_t)
      :abstractmethod:



   .. py:method:: read_variable_surface_air_temperature_at_time(index_t)
      :abstractmethod:



   .. py:method:: read_variable_dew_point_temperature_at_time(index_t)
      :abstractmethod:



   .. py:method:: read_variable_surface_downwards_solar_radiation_at_time(index_t)
      :abstractmethod:



   .. py:method:: read_variable_surface_downwards_thermal_radiation_at_time(index_t)
      :abstractmethod:



   .. py:method:: read_variable_surface_solar_radiation_at_time(index_t)
      :abstractmethod:



   .. py:method:: read_variable_surface_thermal_radiation_at_time(index_t)
      :abstractmethod:



   .. py:method:: read_variable_wind_10m_at_time(index_t)
      :abstractmethod:



