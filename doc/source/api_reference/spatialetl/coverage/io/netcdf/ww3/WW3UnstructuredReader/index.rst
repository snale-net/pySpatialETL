spatialetl.coverage.io.netcdf.ww3.WW3UnstructuredReader
=======================================================

.. py:module:: spatialetl.coverage.io.netcdf.ww3.WW3UnstructuredReader


Classes
-------

.. autoapisummary::

   spatialetl.coverage.io.netcdf.ww3.WW3UnstructuredReader.WW3UnstructuredReader


Module Contents
---------------

.. py:class:: WW3UnstructuredReader(myFile, xsize=None, ysize=None, open_boundaries_mask_to_land=False)

   Bases: :py:obj:`spatialetl.coverage.io.CoverageReader.CoverageReader`


   .. py:attribute:: ncfile


   .. py:attribute:: open_boundaries_mask_to_land
      :value: False



   .. py:method:: close()


   .. py:method:: is_regular_grid()


   .. py:method:: get_x_size()


   .. py:method:: get_y_size()


   .. py:method:: get_t_size()


   .. py:method:: read_axis_t(tmin, tmax, timestamp)


   .. py:method:: read_axis_x(xmin, xmax, ymin, ymax)


   .. py:method:: read_axis_y(xmin, xmax, ymin, ymax)


   .. py:method:: read_variable_2D_sea_binary_mask(xmin, xmax, ymin, ymax)


   .. py:method:: read_variable_bathymetry(xmin, xmax, ymin, ymax)


   .. py:method:: read_variable_bathymetry_at_time(t, xmin, xmax, ymin, ymax)


   .. py:method:: read_variable_sea_surface_height_above_mean_sea_level_at_time(t, xmin, xmax, ymin, ymax)


   .. py:method:: read_variable_sea_surface_wave_significant_height_at_time(t, xmin, xmax, ymin, ymax)


   .. py:method:: read_variable_sea_surface_wave_breaking_height_at_time(t, xmin, xmax, ymin, ymax)


   .. py:method:: read_variable_sea_surface_wave_from_direction_at_time(t, xmin, xmax, ymin, ymax)


   .. py:method:: read_variable_sea_surface_wave_to_direction_at_time(t, xmin, xmax, ymin, ymax)


   .. py:method:: read_variable_sea_surface_wave_mean_period_at_time(t, xmin, xmax, ymin, ymax)


   .. py:method:: read_variable_radiation_pressure_bernouilli_head_at_time(t, xmin, xmax, ymin, ymax)


   .. py:method:: read_variable_sea_surface_wave_energy_flux_to_ocean_at_time(t, xmin, xmax, ymin, ymax)


   .. py:method:: read_variable_sea_surface_wave_peak_frequency_at_time(t, xmin, xmax, ymin, ymax)


   .. py:method:: read_variable_barotropic_sea_water_velocity_at_time(t, xmin, xmax, ymin, ymax)


   .. py:method:: read_variable_atmosphere_momentum_flux_to_waves_at_time(t, xmin, xmax, ymin, ymax)


   .. py:method:: read_variable_waves_momentum_flux_to_ocean_at_time(t, xmin, xmax, ymin, ymax)


   .. py:method:: read_variable_sea_surface_wave_stokes_drift_velocity_at_time(t, xmin, xmax, ymin, ymax)


   .. py:method:: read_variable_wind_10m_at_time(t, xmin, xmax, ymin, ymax)


