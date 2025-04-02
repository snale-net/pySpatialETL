spatialetl.coverage.io.matlab.swan.SWANUnstructuredReader
=========================================================

.. py:module:: spatialetl.coverage.io.matlab.swan.SWANUnstructuredReader


Classes
-------

.. autoapisummary::

   spatialetl.coverage.io.matlab.swan.SWANUnstructuredReader.SWANUnstructuredReader


Module Contents
---------------

.. py:class:: SWANUnstructuredReader(myFile)

   Bases: :py:obj:`spatialetl.coverage.io.CoverageReader.CoverageReader`


   .. py:attribute:: mat


   .. py:attribute:: x_size


   .. py:attribute:: y_size


   .. py:attribute:: times
      :value: []



   .. py:method:: close()


   .. py:method:: is_regular_grid()


   .. py:method:: get_x_size()


   .. py:method:: get_y_size()


   .. py:method:: get_t_size()


   .. py:method:: read_axis_x(xmin, xmax, ymin, ymax)


   .. py:method:: read_axis_y(xmin, xmax, ymin, ymax)


   .. py:method:: read_variable_bathymetry(xmin, xmax, ymin, ymax)


   .. py:method:: read_axis_t(tmin, tmax, timestamp)


   .. py:method:: read_variable_sea_surface_wave_significant_height_at_time(t, xmin, xmax, ymin, ymax)


   .. py:method:: read_variable_sea_surface_wave_mean_period_at_time(t, xmin, xmax, ymin, ymax)


   .. py:method:: read_variable_sea_surface_wave_from_direction_at_time(t, xmin, xmax, ymin, ymax)


   .. py:method:: read_variable_sea_surface_wave_to_direction_at_time(t, xmin, xmax, ymin, ymax)


   .. py:method:: read_variable_wind_10m_at_time(t, xmin, xmax, ymin, ymax)


