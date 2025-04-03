spatialetl.coverage.io.netcdf.mercator.v2020.MERCATORReader
===========================================================

.. py:module:: spatialetl.coverage.io.netcdf.mercator.v2020.MERCATORReader


Classes
-------

.. autoapisummary::

   spatialetl.coverage.io.netcdf.mercator.v2020.MERCATORReader.MERCATORReader


Module Contents
---------------

.. py:class:: MERCATORReader(myFile)

   Bases: :py:obj:`spatialetl.coverage.io.CoverageReader.CoverageReader`


   .. py:method:: close()


   .. py:method:: is_regular_grid()


   .. py:method:: get_x_size()


   .. py:method:: get_y_size()


   .. py:method:: get_z_size()


   .. py:method:: get_t_size()


   .. py:method:: read_axis_x(xmin, xmax, ymin, ymax)


   .. py:method:: read_axis_y(xmin, xmax, ymin, ymax)


   .. py:method:: read_axis_z()


   .. py:method:: read_axis_t(tmin, tmax, timestamp)


   .. py:method:: read_variable_2D_sea_binary_mask(xmin, xmax, ymin, ymax)


   .. py:method:: read_variable_sea_surface_height_above_geoid_at_time(index_t, xmin, xmax, ymin, ymax)


   .. py:method:: read_variable_sea_water_temperature_at_time_and_depth(index_t, index_z, xmin, xmax, ymin, ymax)


   .. py:method:: read_variable_sea_water_salinity_at_time_and_depth(index_t, index_z, xmin, xmax, ymin, ymax)


   .. py:method:: read_variable_baroclinic_sea_water_velocity_at_time_and_depth(index_t, index_z, xmin, xmax, ymin, ymax)


