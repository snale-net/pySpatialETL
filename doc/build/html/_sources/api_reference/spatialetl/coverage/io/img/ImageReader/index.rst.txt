spatialetl.coverage.io.img.ImageReader
======================================

.. py:module:: spatialetl.coverage.io.img.ImageReader


Classes
-------

.. autoapisummary::

   spatialetl.coverage.io.img.ImageReader.ImageReader


Module Contents
---------------

.. py:class:: ImageReader(myFile)

   Bases: :py:obj:`spatialetl.coverage.io.CoverageReader.CoverageReader`


   .. py:attribute:: file


   .. py:attribute:: y_size


   .. py:attribute:: x_size


   .. py:attribute:: y
      :value: None



   .. py:attribute:: x
      :value: None



   .. py:method:: close()


   .. py:method:: is_regular_grid()


   .. py:method:: get_x_size()


   .. py:method:: get_y_size()


   .. py:method:: read_axis_x(xmin, xmax, ymin, ymax)


   .. py:method:: read_axis_y(xmin, xmax, ymin, ymax)


   .. py:method:: read_variable_2D_sea_binary_mask(xmin, xmax, ymin, ymax)


   .. py:method:: read_variable_topography(xmin, xmax, ymin, ymax)


   .. py:method:: read_variable_bathymetry(xmin, xmax, ymin, ymax)


   .. py:method:: pixel2coord(y, x)


