spatialetl.point.io.netcdf.mefoc.MEFOCReader
============================================

.. py:module:: spatialetl.point.io.netcdf.mefoc.MEFOCReader


Classes
-------

.. autoapisummary::

   spatialetl.point.io.netcdf.mefoc.MEFOCReader.MEFOCReader


Module Contents
---------------

.. py:class:: MEFOCReader(myFile, xy, names=None)

   Bases: :py:obj:`spatialetl.point.io.MultiPointReader.MultiPointReader`


   .. py:attribute:: reader


   .. py:attribute:: source_xy_coords
      :value: []



   .. py:attribute:: names
      :value: []



   .. py:attribute:: xy_coords


   .. py:attribute:: xy_values


   .. py:attribute:: meta_data
      :value: ''



   .. py:method:: find_points_coordinates(xy)


   .. py:method:: close()


   .. py:method:: find_point_index(target_lon, target_lat, method='quick', only_mask_value=True)

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



   .. py:method:: is_regular_grid()


   .. py:method:: get_z_size()


   .. py:method:: get_t_size()


   .. py:method:: read_axis_x()


   .. py:method:: read_axis_y()


   .. py:method:: read_axis_t(tmin, tmax, timestamp=0)


   .. py:method:: read_variable_point_names()


   .. py:method:: read_variable_bathymetry()


   .. py:method:: read_variable_sea_surface_height_above_mean_sea_level_at_time(index_t)


   .. py:method:: read_variable_sea_water_column_thickness_at_time(index_t)


   .. py:method:: read_variable_barotropic_sea_water_velocity_at_time(index_t)


   .. py:method:: read_variable_sea_surface_wave_significant_height_at_time(index_t)


   .. py:method:: read_variable_sea_surface_wave_mean_period_at_time(index_t)


   .. py:method:: read_variable_wind_10m_at_time(index_t)


   .. py:method:: read_variable_surface_air_pressure_at_time(index_t)


   .. py:method:: read_variable_rainfall_amount_at_time(index_t)


