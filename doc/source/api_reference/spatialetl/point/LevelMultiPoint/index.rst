spatialetl.point.LevelMultiPoint
================================

.. py:module:: spatialetl.point.LevelMultiPoint


Classes
-------

.. autoapisummary::

   spatialetl.point.LevelMultiPoint.LevelMultiPoint


Module Contents
---------------

.. py:class:: LevelMultiPoint(myReader, zbox=None, resolution_z=None)

   Bases: :py:obj:`spatialetl.point.MultiPoint.MultiPoint`


   .. py:attribute:: DEPTH_DELTA
      :value: 1.0



   .. py:attribute:: VERTICAL_INTERPOLATION_METHOD
      :value: 'nearest'



   .. py:attribute:: vertical_resampling
      :value: False



   .. py:attribute:: source_sigma_coordinate
      :value: False



   .. py:attribute:: target_sigma_coordinate
      :value: False



   .. py:attribute:: source_z_size


   .. py:attribute:: source_axis_z


   .. py:method:: read_axis_z(type='target')


   .. py:method:: get_z_size(type='target')


   .. py:method:: is_sigma_coordinate(type='target')


   .. py:method:: find_level_index(depth, method='fast')

      Retourne l'index de la profondeur la plus proche selon le point le plus proche.
      @type depth : integer ou flottant
      @param depth: Profondeur en mètre souhaitée ou index de la profondeur souhaitée
      @return:  un tableau de l'indice de la couche verticale inférieur la plus proche en chacun point de la grille. z < vert_coord[y,x] et z > vert_coord[y,x]+1.
      Les valeurs masquées valent -999.



   .. py:method:: interpolate_vertical(depth, vert_coord, indexes_z, layers)


   .. py:method:: read_variable_sea_water_temperature_at_depth(depth)


   .. py:method:: read_variable_sea_water_salinity_at_depth(depth)


   .. py:method:: read_variable_sea_water_density_at_depth(depth)


   .. py:method:: read_variable_sea_water_turbidity_at_depth(depth)


   .. py:method:: read_variable_sea_water_electrical_conductivity_at_depth(depth)


