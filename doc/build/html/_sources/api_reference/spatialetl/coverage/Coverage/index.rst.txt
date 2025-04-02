spatialetl.coverage.Coverage
============================

.. py:module:: spatialetl.coverage.Coverage


Classes
-------

.. autoapisummary::

   spatialetl.coverage.Coverage.Coverage


Module Contents
---------------

.. py:class:: Coverage(myReader, bbox=None, resolution_x=None, resolution_y=None)

   Bases: :py:obj:`object`


   La classe Coverage représente une couverture spatiale sur l'horizontale. Les point qui représentent cette couverture
   peuvent être alignés sur une maille régulière (x,y) ou sur une maille non-régulière ((x1,y1),(x2,y2)). En fonction du
   type de maille, les fonctions de lecture des axes retourneront des tableaux à une ou deux dimensions. Pour éviter
   un chargement en mémoire de la totalité du fichier, la coverage contient un pointeur vers un lecteur. Les couches
   sont donc lues à la demande dans le fichier.

   Attention, les axes sont toujours inversés dans les tableaux à cause de NetCDF.
   Soit l'axe y en premier puis l'axe x. Exemple : [y,x]

   @param  myReader: lecteur de fichier


   .. py:attribute:: HORIZONTAL_INTERPOLATION_METHOD
      :value: 'linear'



   .. py:attribute:: HORIZONTAL_OVERLAPING_SIZE
      :value: 2



   .. py:attribute:: reader


   .. py:attribute:: map_mpi
      :value: None



   .. py:attribute:: comm


   .. py:attribute:: size


   .. py:attribute:: rank


   .. py:attribute:: source_regular_grid


   .. py:attribute:: target_regular_grid


   .. py:attribute:: horizontal_resampling
      :value: False



   .. py:attribute:: source_global_x_size


   .. py:attribute:: source_global_y_size


   .. py:attribute:: source_global_axis_x


   .. py:attribute:: source_global_axis_y


   .. py:attribute:: target_global_res_x
      :value: None



   .. py:attribute:: target_global_res_y
      :value: None



   .. py:method:: check_bbox_validity(candidate)


   .. py:method:: check_point_is_inside(target_lon, target_lat, lon, lat, tolerance=5)


   .. py:method:: create_mpi_map()


   .. py:method:: update_mpi_map()


   .. py:method:: read_metadata()

      Lit la metadonnée du fichier si le lecteur contient une fonction read_metadata()


   .. py:method:: get_x_size(type='target', with_overlap=False)


   .. py:method:: get_y_size(type='target', with_overlap=False)


   .. py:method:: is_regular_grid(type='target')

      Retourne vrai si la maille est régulière, sinon faux.
      @return:  vrai si la maille est régulière sinon faux.



   .. py:method:: read_axis_x(type='target', with_overlap=False)

      Retourne les valeurs (souvent la longitude) de l'axe x.
      @return:  un tableau à une ou deux dimensions selon le type de maille des valeurs de l'axe x (souvent la longitude) : [x] ou [y,x].



   .. py:method:: read_axis_y(type='target', with_overlap=False)

      Retourne les valeurs (souvent la latitude) de l'axe y.
      @return:  un tableau à une ou deux dimensions selon le type de maille des valeurs de l'axe y (souvent la latitude) : [x] ou [y,x].



   .. py:method:: find_point_index(target_lon, target_lat, decimal_tolerance=5, method='classic', only_mask_value=True, type='source')

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



   .. py:method:: read_variable_bathymetry()

      Retourne la bathymétrie sur toute la couverture
      @return: un tableau en deux dimensions [y,x].



   .. py:method:: read_variable_topography()

      Retourne la topographie sur toute la couverture
      @return: un tableau en deux dimensions [y,x].



   .. py:method:: read_variable_mesh_size()

      Retourne la taille de la grille sur toute la couverture
      @return: un tableau en deux dimensions [y,x].



   .. py:method:: read_variable_x_mesh_size()

      Retourne la taille de la grille sur l'axe X
      @return: un tableau en deux dimensions [y,x].



   .. py:method:: read_variable_y_mesh_size()

      Retourne la taille de la grille sur l'axe Y
      @return: un tableau en deux dimensions [y,x].



   .. py:method:: read_variable_2D_sea_binary_mask(type='target', with_overlap=False)

      Retourne le masque terre/mer sur toute la couverture
      @return: un tableau en deux dimensions [y,x].
              0 = Terre
              1 = Mer



   .. py:method:: read_variable_Ha()

      Retourne l'amplitude de la réanalyse
      @return: un tableau en deux dimensions [y,x].



