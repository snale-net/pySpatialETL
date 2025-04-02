spatialetl.point.MultiPoint
===========================

.. py:module:: spatialetl.point.MultiPoint


Classes
-------

.. autoapisummary::

   spatialetl.point.MultiPoint.MultiPoint


Functions
---------

.. autoapisummary::

   spatialetl.point.MultiPoint.distance_on_unit_sphere


Module Contents
---------------

.. py:function:: distance_on_unit_sphere(long1, lat1, long2, lat2)

   Calcule la distance en kilomètre entre les deux point. Les coordonnées
   sont données en longitude, latitude.

   @type  long1: number
   @param long1: Coordonnée X du point 1.
   @type  lat1: number
   @param lat1: Coordonnée Y du point 1.
   @type  long2: number
   @param long2: Coordonnée X du point 1.
   @type  lat2: number
   @param lat2: Coordonnée Y du point 1.
   @return:  la ditance en kilomètre en ces deux point sur la Terre.


.. py:class:: MultiPoint(myReader)

   .. py:attribute:: reader


   .. py:attribute:: map_mpi
      :value: None



   .. py:attribute:: comm


   .. py:attribute:: size


   .. py:attribute:: rank


   .. py:attribute:: nb_points


   .. py:attribute:: data_source
      :value: 'Undefined'



   .. py:attribute:: name_station
      :value: 'Undefined'



   .. py:attribute:: x_coord
      :value: 'Undefined'



   .. py:attribute:: y_coord
      :value: 'Undefined'



   .. py:attribute:: vertical_datum
      :value: 'Undefined'



   .. py:attribute:: meta_data
      :value: 'Undefined'



   .. py:method:: read_metadata()

      Lit la metadonnée du fichier si le lecteur contient une fonction read_metadata()


   .. py:method:: get_nb_points()


   .. py:method:: read_axis_x()


   .. py:method:: read_axis_y()


   .. py:method:: find_point_index(target_lon, target_lat, method='classic')

      Retourne le point le plus proche du point donné en paramètre.
      @param target_lon: Coordonnée longitude du point
      @param target_lat: Coordonnée latitude du point
      @param method : Méthode de calcul. "Classic" = On parcourt toute la grille à la recherche du plus prêt.
      @return: un tableau contenant
       [0] : l'index x du point le plus proche
       [1] : la coordonnée en longitude du point le plus proche
       [2] : la coordonnée en latitude point le plus proche
       [3] : la distance du point le plus proche en kilomètre.



   .. py:method:: read_variable_point_names()


   .. py:method:: read_variable_time()

      Read time for all point



   .. py:method:: read_variable_bathymetry()

      Read bathymetry for all point



   .. py:method:: read_variable_bathymetry_at_location(x, y)

      Read bathymetry for a specific point



   .. py:method:: read_variable_sea_surface_temperature()

      Read sea_surface_temperature for all point



   .. py:method:: read_variable_sea_surface_salinity()


   .. py:method:: read_variable_sea_water_pressure_at_sea_water_surface()


   .. py:method:: read_variable_sea_surface_density()


   .. py:method:: read_variable_sea_water_turbidity()


   .. py:method:: read_variable_sea_water_electrical_conductivity()


