spatialetl.point.io.AbstractCoverageReader
==========================================

.. py:module:: spatialetl.point.io.AbstractCoverageReader


Classes
-------

.. autoapisummary::

   spatialetl.point.io.AbstractCoverageReader.AbstractCoverageReader


Module Contents
---------------

.. py:class:: AbstractCoverageReader(myFile, xy, names=None)

   Bases: :py:obj:`spatialetl.point.io.MultiPointReader.MultiPointReader`


   .. py:attribute:: reader
      :value: None



   .. py:attribute:: source_xy_coords
      :value: []



   .. py:attribute:: names
      :value: []



   .. py:attribute:: xy_coords


   .. py:attribute:: xy_values


   .. py:attribute:: meta_data
      :value: ''



   .. py:method:: find_points_coordinates()


   .. py:method:: close()


   .. py:method:: find_point_index(target_lon, target_lat, method='classic', only_mask_value=True)

      Retourne le point le plus proche du point donné en paramètre.
      @param target_lon: Coordonnée longitude du point
      @param target_lat: Coordonnée latitude du point
      @param method : Méthode de calcul. "Classic" = On parcourt toute la grille à la recherche du plus prêt.
      @return: un tableau contenant
       [0] : l'index x du point le plus proche
       [1] : l'index y du point le plus proche
       [2] : la coordonnée en longitude du point le plus proche
       [3] : la coordonnée en latitude point le plus proche
       [4] : la distance du point le plus proche en kilomètre.



   .. py:method:: get_z_size()


   .. py:method:: get_t_size()


   .. py:method:: read_axis_x()


   .. py:method:: read_axis_y()


   .. py:method:: read_axis_z()


   .. py:method:: read_axis_t(tmin, tmax, timestamp=0)


   .. py:method:: read_metadata()


   .. py:method:: read_variable_longitude()


   .. py:method:: read_variable_latitude()


   .. py:method:: read_variable_depth()


   .. py:method:: read_variable_time()


   .. py:method:: read_variable_point_names()


   .. py:method:: read_variable_sea_surface_height_above_mean_sea_level_at_time(index_t)


   .. py:method:: read_variable_sea_water_column_thickness_at_time(index_t)


   .. py:method:: read_variable_sea_surface_temperature_at_time(index_t)


   .. py:method:: read_variable_sea_surface_salinity_at_time(index_t)


   .. py:method:: read_variable_sea_water_velocity_at_sea_water_surface_at_time(index_t)


   .. py:method:: read_variable_sea_water_temperature_at_ground_level_at_time(index_t)


   .. py:method:: read_variable_sea_water_salinity_at_ground_level_at_time(index_t)


   .. py:method:: read_variable_sea_water_velocity_at_ground_level_at_time(index_t)


   .. py:method:: read_variable_bathymetry()


   .. py:method:: read_variable_barotropic_sea_water_velocity_at_time(index_t)


   .. py:method:: read_variable_sea_water_temperature_at_time_and_depth(index_t, index_z)


   .. py:method:: read_variable_sea_water_salinity_at_time_and_depth(index_t, index_z)


   .. py:method:: read_variable_baroclinic_sea_water_velocity_at_time_and_depth(index_t, index_z)


   .. py:method:: read_variable_sea_surface_wave_significant_height_at_time(index_t)


   .. py:method:: read_variable_sea_surface_wave_mean_period_at_time(index_t)


   .. py:method:: read_variable_sea_surface_wave_to_direction_at_time(index_t)


   .. py:method:: read_variable_sea_surface_wave_stokes_drift_velocity_at_time(index_t)


   .. py:method:: read_variable_atmosphere_momentum_flux_to_waves_at_time(index_t)


   .. py:method:: read_variable_waves_momentum_flux_to_ocean_at_time(index_t)


   .. py:method:: read_variable_wind_stress_at_time(index_t)


   .. py:method:: read_variable_wind_10m_at_time(index_t)


