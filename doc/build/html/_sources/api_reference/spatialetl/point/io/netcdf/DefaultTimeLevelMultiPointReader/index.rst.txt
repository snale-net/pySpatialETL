spatialetl.point.io.netcdf.DefaultTimeLevelMultiPointReader
===========================================================

.. py:module:: spatialetl.point.io.netcdf.DefaultTimeLevelMultiPointReader


Classes
-------

.. autoapisummary::

   spatialetl.point.io.netcdf.DefaultTimeLevelMultiPointReader.DefaultTimeLevelMultiPointReader


Module Contents
---------------

.. py:class:: DefaultTimeLevelMultiPointReader(myFile)

   Bases: :py:obj:`spatialetl.point.io.MultiPointReader.MultiPointReader`


   .. py:attribute:: ncfile


   .. py:attribute:: profilmax


   .. py:attribute:: zmax


   .. py:attribute:: tmax


   .. py:method:: get_z_size()


   .. py:method:: read_axis_x()


   .. py:method:: read_axis_y()


   .. py:method:: read_axis_z()


   .. py:method:: read_axis_t(timestamp=0)


   .. py:method:: read_variable_sea_water_temperature_at_time_and_depth(index_t, index_z)


   .. py:method:: read_variable_sea_water_salinity_at_time_and_depth(index_t, index_z)


   .. py:method:: read_variable_sea_water_density_at_time_and_depth(index_t, index_z)


   .. py:method:: read_variable_sea_water_turbidity_at_time_and_depth(index_t, index_z)


   .. py:method:: read_variable_sea_water_electrical_conductivity_at_time_and_depth(index_t, index_z)


