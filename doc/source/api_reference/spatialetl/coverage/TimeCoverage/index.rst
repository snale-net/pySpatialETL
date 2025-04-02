spatialetl.coverage.TimeCoverage
================================

.. py:module:: spatialetl.coverage.TimeCoverage


Classes
-------

.. autoapisummary::

   spatialetl.coverage.TimeCoverage.TimeCoverage


Module Contents
---------------

.. py:class:: TimeCoverage(myReader, bbox=None, resolution_x=None, resolution_y=None, start_time=None, end_time=None, freq=None)

   Bases: :py:obj:`spatialetl.coverage.Coverage.Coverage`


   La classe TimeCoverage est une extension de la classe Coverage.
   Elle rajoute une dimension temporelle à la couverture horizontale classique.



   .. py:attribute:: TIME_DATUM


   .. py:attribute:: TIME_DELTA


   .. py:attribute:: TIME_OVERLAPING_SIZE
      :value: 0



   .. py:attribute:: source_global_t_size


   .. py:attribute:: source_global_axis_t


   .. py:attribute:: temporal_resampling
      :value: False



   .. py:method:: create_mpi_map()


   .. py:method:: update_mpi_map()


   .. py:method:: find_time_index(t, method='fast', domain='source')

      Retourne l'index de la date la plus proche à TIME_DELTA_MIN prêt.
      @type t: datetime ou int
      @param t: date souhaitée ou l'index de la date souhaitée
      @return:  l'index de la date la plus proche à TIME_DELTA_MIN prêt ou une erreur si aucune date n'a pu être trouvée.



   .. py:method:: read_axis_t(type='target', with_overlap=False, timestamp=0)

      Retourne les valeurs de l'axe t.
      @param timestamp: égale 1 si le temps est souhaité en timestamp depuis TIME_DATUM.
      @return:  un tableau à une dimensions [z] au format datetime ou timestamp si timestamp=1.



   .. py:method:: get_t_size(type='target', with_overlap=False)


   .. py:method:: read_variable_2D_sea_binary_mask_at_time(t)

      Retourne le masque à la date souhaitée sur toute la couverture horizontale.
      @type t: datetime ou l'index
      @param t: date souhaitée
      @return: un tableau en deux dimensions [y,x].



   .. py:method:: read_variable_2D_wet_binary_mask_at_time(t)

      Retourne le masque à la date souhaitée sur toute la couverture horizontale.
      @type t: datetime ou l'index
      @param t: date souhaitée
      @return: un tableau en deux dimensions [y,x].



   .. py:method:: read_variable_2D_land_binary_mask_at_time(t)


   .. py:method:: read_variable_sea_surface_height_above_mean_sea_level_at_time(t)


   .. py:method:: read_variable_sea_surface_height_above_geoid_at_time(t)


   .. py:method:: read_variable_sea_water_column_thickness_at_time(t)


   .. py:method:: read_variable_sea_surface_temperature_at_time(t)

      Retourne la temperature de surface à la date souhaitée
      @type t: datetime ou l'index
      @param t: date souhaitée
      @return: un tableau en deux dimensions [y,x].



   .. py:method:: read_variable_sea_surface_salinity_at_time(t)

      Retourne la salinité de surface à la date souhaitée
      @type t: datetime ou l'index
      @param t: date souhaitée
      @return: un tableau en deux dimensions [y,x].



   .. py:method:: read_variable_sea_surface_pressure_at_time(t)

      Retourne la pression à la surface de la mer (sea surface pressure) à la date souhaitée sur toute la couverture horizontale.
      @type t: datetime ou l'index
      @param t: date souhaitée
      @return: un tableau en deux dimensions [y,x].



   .. py:method:: read_variable_sea_surface_density_at_time(t)

      Retourne la densité de l'eau de surface à la date souhaitée
      @type t: datetime ou l'index
      @param t: date souhaitée
      @return: un tableau en deux dimensions [y,x].



   .. py:method:: read_variable_sea_water_turbidity_at_time(t)

      Retourne la turbidité de l'eau de surface à la date souhaitée
      @type t: datetime ou l'index
      @param t: date souhaitée
      @return: un tableau en deux dimensions [y,x].



   .. py:method:: read_variable_sea_water_velocity_at_sea_water_surface_at_time(t)

      Retourne les composantes u,v du courant à la date souhaitée
      @type t: datetime ou l'index
      @param t: date souhaitée
      @return: un tableau en deux dimensions [u_comp,v_comp] contenant chacun deux dimensions [y,x].



   .. py:method:: read_variable_sea_water_temperature_at_ground_level_at_time(t)

      Retourne la temperature de surface à la date souhaitée
      @type t: datetime ou l'index
      @param t: date souhaitée
      @return: un tableau en deux dimensions [y,x].



   .. py:method:: read_variable_sea_water_salinity_at_ground_level_at_time(t)

      Retourne la salinité de surface à la date souhaitée
      @type t: datetime ou l'index
      @param t: date souhaitée
      @return: un tableau en deux dimensions [y,x].



   .. py:method:: read_variable_sea_water_velocity_at_ground_level_at_time(t)

      Retourne les composantes u,v du courant à la date souhaitée
      @type t: datetime ou l'index
      @param t: date souhaitée
      @return: un tableau en deux dimensions [u_comp,v_comp] contenant chacun deux dimensions [y,x].



   .. py:method:: read_variable_barotropic_sea_water_velocity_at_time(t)

      Retourne les composantes u,v du courant à la date souhaitée
      @type t: datetime ou l'index
      @param t: date souhaitée
      @return: un tableau en deux dimensions [u_comp,v_comp] contenant chacun deux dimensions [y,x].



   .. py:method:: read_variable_barotropic_sea_water_speed_at_time(date)


   .. py:method:: read_variable_barotropic_sea_water_from_direction_at_time(date)


   .. py:method:: read_variable_barotropic_sea_water_to_direction_at_time(date)


   .. py:method:: read_variable_sea_surface_wave_significant_height_at_time(t)

      Retourne la hauteur significative des vagues à la date souhaitée sur toute la couverture horizontale.
      @type t: datetime ou l'index
      @param t: date souhaitée
      @return: un tableau en deux dimensions [y,x].



   .. py:method:: read_variable_sea_surface_wave_breaking_height_at_time(t)

      Retourne la hauteur de déferlement des vagues à la date souhaitée sur toute la couverture horizontale.
      @type t: datetime ou l'index
      @param t: date souhaitée
      @return: un tableau en deux dimensions [y,x].



   .. py:method:: read_variable_sea_surface_wave_mean_period_at_time(t)


   .. py:method:: read_variable_sea_surface_wave_peak_period_at_time(t)


   .. py:method:: read_variable_sea_surface_wave_from_direction_at_time(t)


   .. py:method:: read_variable_sea_surface_wave_to_direction_at_time(t)


   .. py:method:: read_variable_sea_surface_wave_stokes_drift_velocity_at_time(t)

      Retourne la dérive de Stokes en surface à la date souhaitée sur toute la couverture horizontale.
      @type t: datetime ou l'index
      @param t: date souhaitée
      @return: un tableau en deux dimensions [y,x].



   .. py:method:: read_variable_radiation_pressure_bernouilli_head_at_time(t)

      Retourne la pression J due aux vagues à la date souhaitée sur toute la couverture horizontale.
      @type t: datetime ou l'index
      @param t: date souhaitée
      @return: un tableau en deux dimensions [y,x].



   .. py:method:: read_variable_sea_surface_wave_energy_flux_to_ocean_at_time(t)

      Retourne la waves_to_ocean_energy_flux à la date souhaitée sur toute la couverture horizontale.
      @type t: datetime ou l'index
      @param t: date souhaitée
      @return: un tableau en deux dimensions [y,x].



   .. py:method:: read_variable_sea_surface_wave_energy_dissipation_at_ground_level_at_time(t)

      Retourne la l'énergie des vagues dissipée par le fond à la date souhaitée sur toute la couverture horizontale.
      @type t: datetime ou l'index
      @param t: date souhaitée
      @return: un tableau en deux dimensions [y,x].



   .. py:method:: read_variable_atmosphere_momentum_flux_to_waves_at_time(t)

      Retourne la composante u du tau atmosphere->vagues à la date souhaitée sur toute la couverture horizontale.
      @type t: datetime ou l'index
      @param t: date souhaitée
      @return: un tableau en deux dimensions [y,x].



   .. py:method:: read_variable_waves_momentum_flux_to_ocean_at_time(t)

      Retourne la composante u du tau vagues->ocean à la date souhaitée sur toute la couverture horizontale.
      @type t: datetime ou l'index
      @param t: date souhaitée
      @return: un tableau en deux dimensions [y,x].



   .. py:method:: read_variable_rainfall_amount_at_time(t)

      Retourne les composantes u,v de rain à la date souhaitée
      @type t: datetime ou l'index
      @param t: date souhaitée
      @return: un tableau en deux dimensions [u_comp,v_comp] contenant chacun deux dimensions [y,x].



   .. py:method:: read_variable_surface_air_pressure_at_time(t)

      Retourne la pression à la surface à la date souhaitée sur toute la couverture horizontale.
      @type t: datetime ou l'index
      @param t: date souhaitée
      @return: un tableau en deux dimensions [y,x].



   .. py:method:: read_variable_sea_surface_air_pressure_at_time(t)

      Retourne la pression à la surface à la date souhaitée sur toute la couverture horizontale.
      @type t: datetime ou l'index
      @param t: date souhaitée
      @return: un tableau en deux dimensions [y,x].



   .. py:method:: read_variable_wind_stress_at_time(t)

      Retourne les composantes u,v de la contrainte de vent à la date souhaitée
      @type t: datetime ou l'index
      @param t: date souhaitée
      @return: un tableau en deux dimensions [u_comp,v_comp] contenant chacun deux dimensions [y,x].



   .. py:method:: read_variable_surface_downward_sensible_heat_flux_at_time(t)

      Retourne les composantes u,v de surface sensible heat flux à la date souhaitée
      @type t: datetime ou l'index
      @param t: date souhaitée
      @return: un tableau en deux dimensions [u_comp,v_comp] contenant chacun deux dimensions [y,x].



   .. py:method:: read_variable_surface_downward_latent_heat_flux_at_time(t)

      Retourne les composantes u,v de surface latente heat flux à la date souhaitée
      @type t: datetime ou l'index
      @param t: date souhaitée
      @return: un tableau en deux dimensions [u_comp,v_comp] contenant chacun deux dimensions [y,x].



   .. py:method:: read_variable_surface_air_temperature_at_time(t)

      Retourne les composantes u,v de surface air temperature à la date souhaitée
      @type t: datetime ou l'index
      @param t: date souhaitée
      @return: un tableau en deux dimensions [u_comp,v_comp] contenant chacun deux dimensions [y,x].



   .. py:method:: read_variable_dew_point_temperature_at_time(t)

      Retourne les composantes u,v de dewpoint temperature à la date souhaitée
      @type t: datetime ou l'index
      @param t: date souhaitée
      @return: un tableau en deux dimensions [u_comp,v_comp] contenant chacun deux dimensions [y,x].



   .. py:method:: read_variable_surface_downward_solar_radiation_at_time(t)

      Retourne les composantes u,v de surface solar radiation downwards à la date souhaitée
      @type t: datetime ou l'index
      @param t: date souhaitée
      @return: un tableau en deux dimensions [u_comp,v_comp] contenant chacun deux dimensions [y,x].



   .. py:method:: read_variable_surface_downward_thermal_radiation_at_time(t)

      Retourne les composantes u,v de surface thermal radiation downwards à la date souhaitée
      @type t: datetime ou l'index
      @param t: date souhaitée
      @return: un tableau en deux dimensions [u_comp,v_comp] contenant chacun deux dimensions [y,x].



   .. py:method:: read_variable_surface_solar_radiation_at_time(t)

      Retourne les composantes u,v de surface solar radiation à la date souhaitée
      @type t: datetime ou l'index
      @param t: date souhaitée
      @return: un tableau en deux dimensions [u_comp,v_comp] contenant chacun deux dimensions [y,x].



   .. py:method:: read_variable_surface_thermal_radiation_at_time(t)

      Retourne les composantes u,v de surface thermal radiation à la date souhaitée
      @type t: datetime ou l'index
      @param t: date souhaitée
      @return: un tableau en deux dimensions [u_comp,v_comp] contenant chacun deux dimensions [y,x].



   .. py:method:: read_variable_wind_10m_at_time(t)

      Retourne les composantes u,v du vent à la date souhaitée
      @type t: datetime ou l'index
      @param t: date souhaitée
      @return: un tableau en deux dimensions [u_comp,v_comp] contenant chacun deux dimensions [y,x].



   .. py:method:: read_variable_wind_speed_10m_at_time(date)


   .. py:method:: read_variable_wind_from_direction_10m_at_time(date)


   .. py:method:: read_variable_wind_to_direction_10m_at_time(date)


