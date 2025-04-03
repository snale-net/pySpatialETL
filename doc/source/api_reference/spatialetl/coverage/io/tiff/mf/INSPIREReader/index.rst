spatialetl.coverage.io.tiff.mf.INSPIREReader
============================================

.. py:module:: spatialetl.coverage.io.tiff.mf.INSPIREReader


Classes
-------

.. autoapisummary::

   spatialetl.coverage.io.tiff.mf.INSPIREReader.INSPIREReader


Module Contents
---------------

.. py:class:: INSPIREReader(myFile)

   Bases: :py:obj:`spatialetl.coverage.io.CoverageReader.CoverageReader`


   .. py:attribute:: VARIABLES
      :value: ['BRIGHTNESS_TEMPERATURE__GROUND_OR_WATER_SURFACE',...



   .. py:attribute:: VARIABLES_PT3H
      :value: ['TOTAL_WATER_PRECIPITATION__GROUND_OR_WATER_SURFACE',...



   .. py:attribute:: times
      :value: []



   .. py:attribute:: variable_files


   .. py:attribute:: last_opened_t_index
      :value: 0



   .. py:attribute:: tifffile


   .. py:attribute:: y
      :value: None



   .. py:attribute:: x
      :value: None



   .. py:method:: open_file(varname, index_t)


   .. py:method:: close()


   .. py:method:: is_regular_grid()


   .. py:method:: pixel2coord(y, x)


   .. py:method:: get_x_size()


   .. py:method:: get_y_size()


   .. py:method:: get_t_size()


   .. py:method:: read_axis_t(tmin, tmax, timestamp)


   .. py:method:: read_axis_x(xmin, xmax, ymin, ymax)


   .. py:method:: read_axis_y(xmin, xmax, ymin, ymax)


   .. py:method:: read_variable_2D_sea_binary_mask(xmin, xmax, ymin, ymax)


   .. py:method:: read_variable_2D_sea_binary_mask_at_time(intex_t, xmin, xmax, ymin, ymax)


   .. py:method:: read_variable_2D_land_binary_mask(xmin, xmax, ymin, ymax)


   .. py:method:: read_variable_2D_land_binary_mask_at_time(intex_t, xmin, xmax, ymin, ymax)


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


