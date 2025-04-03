spatialetl.coverage.io.netcdf.mefoc.MEFOCReader
===============================================

.. py:module:: spatialetl.coverage.io.netcdf.mefoc.MEFOCReader


Classes
-------

.. autoapisummary::

   spatialetl.coverage.io.netcdf.mefoc.MEFOCReader.MEFOCReader


Module Contents
---------------

.. py:class:: MEFOCReader(myFile)

   Bases: :py:obj:`spatialetl.coverage.io.CoverageReader.CoverageReader`


   .. py:method:: close()


   .. py:method:: is_regular_grid()


   .. py:method:: get_x_size()


   .. py:method:: get_y_size()


   .. py:method:: get_z_size()


   .. py:method:: get_t_size()


   .. py:method:: read_axis_t(tmin, tmax, timestamp=0)


   .. py:method:: read_axis_x(xmin, xmax, ymin, ymax)


   .. py:method:: read_axis_y(xmin, xmax, ymin, ymax)


   .. py:method:: read_variable_bathymetry(xmin, xmax, ymin, ymax)


   .. py:method:: read_variable_sea_surface_height_above_mean_sea_level_at_time(index_t, xmin, xmax, ymin, ymax)


   .. py:method:: read_variable_sea_water_column_thickness_at_time(index_t, xmin, xmax, ymin, ymax)


   .. py:method:: read_variable_barotropic_sea_water_velocity_at_time(index_t, xmin, xmax, ymin, ymax)


   .. py:method:: read_variable_sea_surface_wave_significant_height_at_time(index_t, xmin, xmax, ymin, ymax)


   .. py:method:: read_variable_sea_surface_wave_mean_period_at_time(index_t, xmin, xmax, ymin, ymax)


   .. py:method:: read_variable_sea_surface_wave_to_direction_at_time(index_t, xmin, xmax, ymin, ymax)


   .. py:method:: read_variable_eastward_wind_10m(xmin, xmax, ymin, ymax)


   .. py:method:: read_variable_northward_wind_10m(xmin, xmax, ymin, ymax)


   .. py:method:: read_variable_wind_10m_at_time(index_t, xmin, xmax, ymin, ymax)


