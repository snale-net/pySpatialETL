spatialetl.coverage.io.netcdf.symphonie.old.SYMPHONIEReader
===========================================================

.. py:module:: spatialetl.coverage.io.netcdf.symphonie.old.SYMPHONIEReader


Classes
-------

.. autoapisummary::

   spatialetl.coverage.io.netcdf.symphonie.old.SYMPHONIEReader.SYMPHONIEReader


Module Contents
---------------

.. py:class:: SYMPHONIEReader(myGrid, myFile)

   Bases: :py:obj:`spatialetl.coverage.io.netcdf.symphonie.SYMPHONIEReader.SYMPHONIEReader`


   La classe SymphonieReader permet de lire les données du format Symphonie

   @param  myGrid : lien vers le fichier de grille (que l'on trouve dans le RDIR/tmp/grid.nc)
   @param myFile : lien vers le fichier de données (que l'on trouve dans GRAPHIQUES)


   .. py:method:: read_variable_current_at_time_and_depth(index_t, index_z)


