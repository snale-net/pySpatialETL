spatialetl.coverage.TimeLevelCoverage
=====================================

.. py:module:: spatialetl.coverage.TimeLevelCoverage


Classes
-------

.. autoapisummary::

   spatialetl.coverage.TimeLevelCoverage.TimeLevelCoverage


Module Contents
---------------

.. py:class:: TimeLevelCoverage(myReader, bbox=None, resolution_x=None, resolution_y=None, zbox=None, resolution_z=None, start_time=None, end_time=None, freq=None)

   Bases: :py:obj:`spatialetl.coverage.LevelCoverage.LevelCoverage`, :py:obj:`spatialetl.coverage.TimeCoverage.TimeCoverage`


   La classe TimeLevelCoverage est une extension de la classe Coverage, LevelCoverage, TimeCoverage.
   Elle rajoute les dimensions temporelle et verticale à la couverture horizontale classique.



   .. py:attribute:: data_temp


   .. py:attribute:: layers_temp


   .. py:method:: read_variable_sea_water_temperature_at_time_and_depth(time, depth)

      Retourne la salinité à la date souhaitée et au niveau souhaité sur toute la couverture horizontale.
      @type time: datetime ou l'index
      @param time: date souhaitée
      @type depth: profondeur en mètre (float) ou index (integer)
      @param depth: profondeur souhaitée. Si le z est un entier, on considère qu'il s'agit de l'index,
      si c'est un flottant on considère qu'il s'agit d'une profondeur
      @return: un tableau en deux dimensions [y,x].



   .. py:method:: read_variable_sea_water_salinity_at_time_and_depth(time, depth)

      Retourne la salinité à la date souhaitée et au niveau souhaité sur toute la couverture horizontale.
      @type time: datetime ou l'index
      @param time: date souhaitée
      @type depth: profondeur en mètre (float) ou index (integer)
      @param depth: profondeur souhaitée. Si le z est un entier, on considère qu'il s'agit de l'index,
      si c'est un flottant on considère qu'il s'agit d'une profondeur
      @return: un tableau en deux dimensions [y,x].



   .. py:method:: read_variable_baroclinic_sea_water_velocity_at_time_and_depth(time, depth)

      Retourne les composantes u,v du courant à la date souhaitée et au niveau souhaité sur toute la couverture horizontale.
      @type time: datetime ou l'index
      @param time: date souhaitée
      @type depth: profondeur en mètre (float) ou index (integer)
      @param depth: profondeur souhaitée. Si le z est un entier, on considère qu'il s'agit de l'index,
      si c'est un flottant on considère qu'il s'agit d'une profondeur
      @return: un tableau en deux dimensions [u_comp,v_comp] contenant chacun deux dimensions [y,x].



