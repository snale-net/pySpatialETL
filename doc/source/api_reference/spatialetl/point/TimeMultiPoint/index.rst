spatialetl.point.TimeMultiPoint
===============================

.. py:module:: spatialetl.point.TimeMultiPoint


Classes
-------

.. autoapisummary::

   spatialetl.point.TimeMultiPoint.TimeMultiPoint


Module Contents
---------------

.. py:class:: TimeMultiPoint(myReader, start_time=None, end_time=None, freq=None, time_range=None)

   Bases: :py:obj:`spatialetl.point.MultiPoint.MultiPoint`


   .. py:attribute:: TIME_DATUM


   .. py:attribute:: TIME_DELTA


   .. py:attribute:: TIME_INTERPOLATION_METHOD
      :value: 'nearest'



   .. py:attribute:: TIME_OVERLAPING_SIZE
      :value: 2



   .. py:attribute:: source_global_t_size


   .. py:attribute:: source_global_axis_t


   .. py:attribute:: temporal_resampling
      :value: False



   .. py:method:: create_mpi_map()


   .. py:method:: update_mpi_map()


   .. py:method:: read_axis_t(type='target', with_overlap=False, timestamp=0)

      Retourne les valeurs de l'axe t.
      @param timestamp: égale 1 si le temps est souhaité en timestamp depuis TIME_DATUM.
      @return:  un tableau à une dimensions [z] au format datetime ou timestamp si timestamp=1.



   .. py:method:: get_t_size(type='target', with_overlap=False)


   .. py:method:: find_time_index(t)

      Retourne l'index de la date la plus proche à TIME_DELTA_MIN prêt.
      @type t: datetime ou int
      @param t: date souhaitée ou l'index de la date souhaitée
      @return:  l'index de la date la plus proche à TIME_DELTA_MIN prêt ou une erreur si aucune date n'a pu être trouvée.



   .. py:method:: interpolate_all_times(values)


   .. py:method:: interpolate_time(date, indexes_t, layers)


   .. py:method:: read_variable_longitude_at_time(date)


   .. py:method:: read_variable_latitude_at_time(date)


   .. py:method:: read_variable_sea_surface_height_above_mean_sea_level_at_time(date)


   .. py:method:: read_variable_sea_surface_height_above_geoid_at_time(date)


   .. py:method:: read_variable_sea_water_column_thickness_at_time(date)


   .. py:method:: read_variable_sea_surface_density_at_time(date)


   .. py:method:: read_variable_sea_water_turbidity_at_time(date)


   .. py:method:: read_variable_sea_water_electrical_conductivity_at_time(date)


   .. py:method:: read_variable_barotropic_sea_water_velocity_at_time(date)


   .. py:method:: read_variable_barotropic_sea_water_speed_at_time(date)


   .. py:method:: read_variable_barotropic_sea_water_from_direction_at_time(date)


   .. py:method:: read_variable_barotropic_sea_water_to_direction_at_time(date)


   .. py:method:: read_variable_sea_surface_temperature_at_time(date)


   .. py:method:: read_variable_sea_surface_salinity_at_time(date)


   .. py:method:: read_variable_sea_water_pressure_at_sea_water_surface_at_time(date)


   .. py:method:: read_variable_sea_water_velocity_at_sea_water_surface_at_time(date)


   .. py:method:: read_variable_sea_water_speed_at_sea_water_surface_at_time(date)


   .. py:method:: read_variable_sea_water_from_direction_at_sea_water_surface_at_time(date)


   .. py:method:: read_variable_sea_water_to_direction_at_sea_water_surface_at_time(date)


   .. py:method:: read_variable_sea_water_temperature_at_ground_level_at_time(date)


   .. py:method:: read_variable_sea_water_salinity_at_ground_level_at_time(date)


   .. py:method:: read_variable_sea_water_velocity_at_ground_level_at_time(date)


   .. py:method:: read_variable_sea_water_speed_at_ground_level_at_time(date)


   .. py:method:: read_variable_sea_water_from_direction_at_ground_level_at_time(date)


   .. py:method:: read_variable_sea_water_to_direction_at_ground_level_at_time(date)


   .. py:method:: read_variable_sea_surface_wave_significant_height_at_time(date)


   .. py:method:: read_variable_sea_surface_wave_mean_period_at_time(date)


   .. py:method:: read_variable_sea_surface_wave_to_direction_at_time(date)


   .. py:method:: read_variable_water_volume_transport_into_sea_water_from_rivers_at_time(date)


   .. py:method:: read_variable_surface_air_pressure_at_time(date)


   .. py:method:: read_variable_sea_surface_air_pressure_at_time(date)


   .. py:method:: read_variable_surface_downward_sensible_heat_flux_at_time(date)


   .. py:method:: read_variable_rainfall_amount_at_time(date)


   .. py:method:: read_variable_wind_stress_at_time(date)


   .. py:method:: read_variable_wind_stress_stress_at_time(date)


   .. py:method:: read_variable_wind_stress_from_direction_at_time(date)


   .. py:method:: read_variable_wind_stress_to_direction_at_time(date)


   .. py:method:: read_variable_surface_downward_latent_heat_flux_at_time(date)


   .. py:method:: read_variable_surface_air_temperature_at_time(date)


   .. py:method:: read_variable_dew_point_temperature_at_time(date)


   .. py:method:: read_variable_surface_downward_solar_radiation_at_time(date)


   .. py:method:: read_variable_surface_downward_thermal_radiation_at_time(date)


   .. py:method:: read_variable_surface_solar_radiation_at_time(date)


   .. py:method:: read_variable_surface_thermal_radiation_at_time(date)


   .. py:method:: read_variable_wind_10m_at_time(date)


   .. py:method:: read_variable_wind_speed_10m_at_time(date)


   .. py:method:: read_variable_wind_from_direction_10m_at_time(date)


   .. py:method:: read_variable_wind_to_direction_10m_at_time(date)


   .. py:method:: read_variable_sea_surface_temperature()

      Read sea_surface_temperature for all point



   .. py:method:: read_variable_sea_surface_salinity()


   .. py:method:: read_variable_sea_water_pressure_at_sea_water_surface()


   .. py:method:: read_variable_sea_surface_density()


   .. py:method:: read_variable_sea_water_turbidity()


   .. py:method:: read_variable_sea_water_electrical_conductivity()


