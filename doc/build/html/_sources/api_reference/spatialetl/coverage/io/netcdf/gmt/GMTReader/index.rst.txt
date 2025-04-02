spatialetl.coverage.io.netcdf.gmt.GMTReader
===========================================

.. py:module:: spatialetl.coverage.io.netcdf.gmt.GMTReader


Classes
-------

.. autoapisummary::

   spatialetl.coverage.io.netcdf.gmt.GMTReader.GMTReader


Module Contents
---------------

.. py:class:: GMTReader(myFile)

   Bases: :py:obj:`spatialetl.coverage.io.CoverageReader.CoverageReader`


   .. py:attribute:: ncfile


   .. py:method:: close()


   .. py:method:: is_regular_grid()


   .. py:method:: get_x_size()


   .. py:method:: get_y_size()


   .. py:method:: read_axis_x(xmin, xmax, ymin, ymax)


   .. py:method:: read_axis_y(xmin, xmax, ymin, ymax)


   .. py:method:: read_variable_longitude(xmin, xmax, ymin, ymax)


   .. py:method:: read_variable_latitude(xmin, xmax, ymin, ymax)


   .. py:method:: read_variable_2D_sea_binary_mask(xmin, xmax, ymin, ymax)


   .. py:method:: read_variable_bathymetry(xmin, xmax, ymin, ymax)


