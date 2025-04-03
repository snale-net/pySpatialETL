spatialetl.point.io.ascii.DefaultTimePointReader
================================================

.. py:module:: spatialetl.point.io.ascii.DefaultTimePointReader


Classes
-------

.. autoapisummary::

   spatialetl.point.io.ascii.DefaultTimePointReader.DefaultTimePointReader


Module Contents
---------------

.. py:class:: DefaultTimePointReader(myFilename, names, colsNumber, varNames, checkOverlapping=False)

   Bases: :py:obj:`spatialetl.point.io.MultiPointReader.MultiPointReader`


   .. py:attribute:: names


   .. py:attribute:: header
      :value: 0



   .. py:attribute:: x
      :value: ['Undefinied']



   .. py:attribute:: y
      :value: ['Undefinied']



   .. py:method:: read_axis_x()


   .. py:method:: read_axis_y()


   .. py:method:: read_axis_t(timestamp=0)


   .. py:method:: read_metadata()


   .. py:method:: read_variable_wind_10m_at_time(index_t)


   .. py:method:: read_variable_surface_air_pressure_at_time(index_t)


   .. py:method:: read_variable_rainfall_amount_at_time(index_t)


   .. py:method:: read_variable_point_names()


   .. py:method:: read_variable_water_volume_transport_into_sea_water_from_rivers_at_time(index_t)


   .. py:method:: read_variable_sea_water_pressure_at_sea_water_surface_at_time(index_t)


   .. py:method:: read_variable_sea_surface_height_above_mean_sea_level_at_time(index_t)


