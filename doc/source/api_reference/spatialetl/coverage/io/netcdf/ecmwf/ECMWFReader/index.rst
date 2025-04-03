spatialetl.coverage.io.netcdf.ecmwf.ECMWFReader
===============================================

.. py:module:: spatialetl.coverage.io.netcdf.ecmwf.ECMWFReader


Classes
-------

.. autoapisummary::

   spatialetl.coverage.io.netcdf.ecmwf.ECMWFReader.ECMWFReader


Module Contents
---------------

.. py:class:: ECMWFReader(myFile)

   Bases: :py:obj:`spatialetl.coverage.io.CoverageReader.CoverageReader`


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


   .. py:method:: read_variable_2D_land_binary_mask(xmin, xmax, ymin, ymax)


   .. py:method:: read_variable_2D_sea_binary_mask(xmin, xmax, ymin, ymax)


   .. py:method:: read_variable_3D_sea_binary_mask_at_time(index_t, xmin, xmax, ymin, ymax)


   .. py:method:: read_variable_3D_land_binary_mask_at_time(index_t, xmin, xmax, ymin, ymax)


   .. py:method:: read_variable_rainfall_amount_at_time(index_t, xmin, xmax, ymin, ymax)


   .. py:method:: read_variable_surface_air_pressure_at_time(index_t, xmin, xmax, ymin, ymax)


   .. py:method:: read_variable_sea_surface_air_pressure_at_time(index_t, xmin, xmax, ymin, ymax)


   .. py:method:: read_variable_surface_downward_sensible_heat_flux_at_time(index_t, xmin, xmax, ymin, ymax)


   .. py:method:: read_variable_surface_downward_latent_heat_flux_at_time(index_t, xmin, xmax, ymin, ymax)


   .. py:method:: read_variable_surface_air_temperature_at_time(index_t, xmin, xmax, ymin, ymax)


   .. py:method:: read_variable_dew_point_temperature_at_time(index_t, xmin, xmax, ymin, ymax)


   .. py:method:: read_variable_surface_downward_solar_radiation_at_time(index_t, xmin, xmax, ymin, ymax)


   .. py:method:: read_variable_surface_downward_thermal_radiation_at_time(index_t, xmin, xmax, ymin, ymax)


   .. py:method:: read_variable_surface_solar_radiation_at_time(index_t, xmin, xmax, ymin, ymax)


   .. py:method:: read_variable_surface_thermal_radiation_at_time(index_t, xmin, xmax, ymin, ymax)


   .. py:method:: read_variable_wind_10m_at_time(index_t, xmin, xmax, ymin, ymax)


