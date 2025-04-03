spatialetl.coverage.io.netcdf.ecmwf.AcademicECMWFWriter
=======================================================

.. py:module:: spatialetl.coverage.io.netcdf.ecmwf.AcademicECMWFWriter


Classes
-------

.. autoapisummary::

   spatialetl.coverage.io.netcdf.ecmwf.AcademicECMWFWriter.AcademicECMWFWriter


Module Contents
---------------

.. py:class:: AcademicECMWFWriter(myFile, lon, lat, times, wind_speed=0, wind_from_direction_angle=0, surface_air_pressure=1013.25, sea_surface_air_pressure=1013.25, surface_air_temperature=283.15, dewpoint_temperature=283.15, surface_downward_sensible_heat_flux=0, surface_downward_latent_heat_flux=0, surface_downward_solar_radiation=0, surface_downward_thermal_radiation=0, surface_solar_radiation=0, surface_thermal_radiation=0, total_rain=0, update=False)

   Bases: :py:obj:`spatialetl.coverage.io.CoverageWriter.CoverageWriter`


   .. py:attribute:: x_axis


   .. py:attribute:: y_axis


   .. py:attribute:: t_axis


   .. py:attribute:: wind_speed
      :value: 0



   .. py:attribute:: wind_from_direction_angle
      :value: 0



   .. py:attribute:: surface_air_pressure
      :value: 1013.25



   .. py:attribute:: sea_surface_air_pressure
      :value: 1013.25



   .. py:attribute:: surface_downward_sensible_heat_flux
      :value: 0



   .. py:attribute:: surface_downward_latent_heat_flux
      :value: 0



   .. py:attribute:: surface_air_temperature
      :value: 283.15



   .. py:attribute:: dewpoint_temperature
      :value: 283.15



   .. py:attribute:: surface_downward_solar_radiation
      :value: 0



   .. py:attribute:: surface_solar_radiation
      :value: 0



   .. py:attribute:: surface_thermal_radiation
      :value: 0



   .. py:attribute:: surface_downward_thermal_radiation
      :value: 0



   .. py:attribute:: total_rain
      :value: 0



   .. py:attribute:: ncfile
      :value: None



   .. py:attribute:: update
      :value: False



   .. py:method:: close()


   .. py:method:: write_variable_3D_mask()


   .. py:method:: write_variable_surface_pressure()


   .. py:method:: write_variable_sea_surface_air_pressure()


   .. py:method:: write_variable_wind()


   .. py:method:: write_variable_surface_downward_sensible_heat_flux()


   .. py:method:: write_variable_surface_downward_latent_heat_flux()


   .. py:method:: write_variable_surface_air_temperature()


   .. py:method:: write_variable_dewpoint_temperature()


   .. py:method:: write_variable_surface_downward_solar_radiation()


   .. py:method:: write_variable_surface_downward_thermal_radiation()


   .. py:method:: write_variable_surface_solar_radiation()


   .. py:method:: write_variable_surface_thermal_radiation()


   .. py:method:: write_variable_rainfall_amount()


