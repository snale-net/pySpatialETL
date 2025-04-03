spatialetl.coverage.io.grib.meteofrance.MFWAMReader
===================================================

.. py:module:: spatialetl.coverage.io.grib.meteofrance.MFWAMReader


Classes
-------

.. autoapisummary::

   spatialetl.coverage.io.grib.meteofrance.MFWAMReader.MFWAMReader


Module Contents
---------------

.. py:class:: MFWAMReader(myFile)

   Bases: :py:obj:`spatialetl.coverage.io.CoverageReader.CoverageReader`


   .. py:attribute:: t_size
      :value: 1



   .. py:attribute:: times


   .. py:attribute:: new_lon


   .. py:attribute:: new_indexes


   .. py:method:: find_time_and_step(index_t)


   .. py:method:: close()


   .. py:method:: is_regular_grid()


   .. py:method:: get_x_size()


   .. py:method:: get_y_size()


   .. py:method:: get_t_size()


   .. py:method:: read_axis_x(xmin, xmax, ymin, ymax)


   .. py:method:: read_axis_y(xmin, xmax, ymin, ymax)


   .. py:method:: read_axis_t(tmin, tmax, timestamp)


   .. py:method:: read_variable_longitude(xmin, xmax, ymin, ymax)


   .. py:method:: read_variable_latitude(xmin, xmax, ymin, ymax)


   .. py:method:: read_variable_time(tmin, tmax)


   .. py:method:: read_variable_sea_surface_wave_significant_height_at_time(index_t, xmin, xmax, ymin, ymax)


