spatialetl.coverage.io.netcdf.mercator.v2015.MERCATORReader
===========================================================

.. py:module:: spatialetl.coverage.io.netcdf.mercator.v2015.MERCATORReader


Classes
-------

.. autoapisummary::

   spatialetl.coverage.io.netcdf.mercator.v2015.MERCATORReader.MERCATORReader


Module Contents
---------------

.. py:class:: MERCATORReader(m, d, t, u, v)

   Bases: :py:obj:`spatialetl.coverage.io.CoverageReader.CoverageReader`


   .. py:attribute:: mask


   .. py:attribute:: grid2D


   .. py:method:: close()


   .. py:method:: is_regular_grid()


   .. py:method:: get_x_size()


   .. py:method:: get_y_size()


   .. py:method:: get_z_size()


   .. py:method:: get_t_size()


   .. py:method:: read_axis_t(tmin, tmax, timestamp)

      Attention si gridT, U,V,2D ont un time_counter different



   .. py:method:: read_axis_x(xmin, xmax, ymin, ymax)


   .. py:method:: read_axis_y(xmin, xmax, ymin, ymax)


   .. py:method:: read_axis_z()


   .. py:method:: read_variable_2D_sea_binary_mask(xmin, xmax, ymin, ymax)


   .. py:method:: read_variable_3D_sea_binary_mask(xmin, xmax, ymin, ymax)


   .. py:method:: read_variable_4D_sea_binary_mask(xmin, xmax, ymin, ymax)


   .. py:method:: read_variable_sea_surface_height_above_sea_level_at_time(t, xmin, xmax, ymin, ymax)


   .. py:method:: read_variable_baroclinic_sea_water_velocity_at_time_and_depth(index_t, index_z, xmin, xmax, ymin, ymax)


