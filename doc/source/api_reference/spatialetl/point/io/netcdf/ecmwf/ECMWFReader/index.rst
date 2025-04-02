spatialetl.point.io.netcdf.ecmwf.ECMWFReader
============================================

.. py:module:: spatialetl.point.io.netcdf.ecmwf.ECMWFReader


Classes
-------

.. autoapisummary::

   spatialetl.point.io.netcdf.ecmwf.ECMWFReader.ECMWFReader


Module Contents
---------------

.. py:class:: ECMWFReader(myFile, xy, names=None)

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



   .. py:method:: get_t_size()


   .. py:method:: read_axis_x()


   .. py:method:: read_axis_y()


   .. py:method:: read_axis_t(tmin, tmax, timestamp=0)


   .. py:method:: read_metadata()


   .. py:method:: read_variable_point_names()


   .. py:method:: read_variable_wind_10m_at_time(index_t)


   .. py:method:: read_variable_surface_air_pressure_at_time(index_t)


   .. py:method:: read_variable_sea_surface_air_pressure_at_time(index_t)


   .. py:method:: read_variable_rainfall_amount_at_time(index_t)


   .. py:method:: read_variable_surface_downward_sensible_heat_flux_at_time(index_t)


   .. py:method:: read_variable_surface_downward_latent_heat_flux_at_time(index_t)


   .. py:method:: read_variable_surface_downward_solar_radiation_at_time(index_t)


   .. py:method:: read_variable_surface_downward_thermal_radiation_at_time(index_t)


   .. py:method:: read_variable_surface_solar_radiation_at_time(index_t)


   .. py:method:: read_variable_surface_thermal_radiation_at_time(index_t)


   .. py:method:: read_variable_surface_air_temperature_at_time(index_t)


   .. py:method:: read_variable_dew_point_temperature_at_time(index_t)


