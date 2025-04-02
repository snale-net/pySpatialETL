spatialetl.coverage.LevelCoverage
=================================

.. py:module:: spatialetl.coverage.LevelCoverage


Classes
-------

.. autoapisummary::

   spatialetl.coverage.LevelCoverage.LevelCoverage


Module Contents
---------------

.. py:class:: LevelCoverage(myReader, bbox=None, resolution_x=None, resolution_y=None, zbox=None, resolution_z=None)

   Bases: :py:obj:`spatialetl.coverage.Coverage.Coverage`


   La classe LevelCoverage est une extension de la classe Coverage.
   Elle rajoute une dimension verticale à la couverture horizontale classique.


   .. py:attribute:: DEPTH_DELTA
      :value: 1.0



   .. py:attribute:: VERTICAL_INTERPOLATION_METHOD
      :value: 'linear'



   .. py:attribute:: vertical_resampling
      :value: False



   .. py:attribute:: source_sigma_coordinate
      :value: False



   .. py:attribute:: target_sigma_coordinate
      :value: False



   .. py:attribute:: source_global_z_size


   .. py:attribute:: source_global_axis_z


   .. py:attribute:: depth_weight


   .. py:method:: read_axis_z(type='target', with_horizontal_overlap=False)

      Retourne les valeurs (souvent en mètre) de l'axe z.
      @return:  si la grille est en coordonnée sigma alors un tableau à trois dimensions [z,y,x] est retourné sinon
      un tableau une dimension [z].



   .. py:method:: is_sigma_coordinate(type='target')

      Retourne vrai si la grille verticale est en coordonnée sigma, sinon faux.
      @return:  vrai si la grille verticale est en coordonnée sigma sinon faux.



   .. py:method:: get_z_size(type='target')

      Retourne la taille de l'axe z.
      @return:  un entier correspondant à la taille de l'axe z.



   .. py:method:: find_level_index(depth, method='fast')

      Retourne l'index de la profondeur la plus proche selon le point le plus proche.
      @type depth : integer ou flottant
      @param depth: Profondeur en mètre souhaitée ou index de la profondeur souhaitée
      @return:  un tableau de l'indice de la couche verticale inférieur la plus proche en chacun point de la grille. z < vert_coord[y,x] et z > vert_coord[y,x]+1.
      Les valeurs masquées valent -999.



   .. py:method:: read_variable_3D_sea_binary_mask()

      Retourne le masque terre/mer sur toute la couverture selon la profondeur z
      @return: un tableau en deux dimensions [z,y,x].
              0 = Terre
              1 = Mer



   .. py:method:: read_variable_depth_at_depth(depth)


