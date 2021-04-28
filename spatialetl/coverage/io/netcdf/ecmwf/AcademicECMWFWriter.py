#! /usr/bin/env python2.7
# -*- coding: utf-8 -*-
#
# pySpatialETL is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# any later version.
#
# pySpatialETL is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# Author : Fabien Rétif - fabien.retif@zoho.com
#
import math

from netCDF4 import Dataset
from netCDF4 import date2num
from numpy import float32
from numpy import float64

from spatialetl.coverage.io.CoverageWriter import CoverageWriter


class AcademicECMWFWriter (CoverageWriter):

    def __init__(self, myFile,
                 lon,
                 lat,
                 times,
                 wind_speed=0,
                 wind_from_direction_angle=0,
                 surface_air_pressure=1013.25,  # Pression MSL
                 sea_surface_air_pressure=1013.25, # Pression MSL
                 surface_air_temperature=283.15, # 10°C
                 dewpoint_temperature=283.15, # 10°C
                 surface_downward_sensible_heat_flux=0,
                 surface_downward_latent_heat_flux=0,
                 surface_downward_solar_radiation=0,
                 surface_downward_thermal_radiation=0,
                 surface_solar_radiation=0,
                 surface_thermal_radiation=0,
                 total_rain=0,
                 update=False):
        CoverageWriter.__init__(self,None,myFile);
        self.x_axis = lon;
        self.y_axis= lat;
        self.t_axis = times
        self.wind_speed = wind_speed
        self.wind_from_direction_angle = wind_from_direction_angle
        self.surface_air_pressure = surface_air_pressure
        self.sea_surface_air_pressure = sea_surface_air_pressure
        self.surface_downward_sensible_heat_flux = surface_downward_sensible_heat_flux
        self.surface_downward_latent_heat_flux = surface_downward_latent_heat_flux
        self.surface_air_temperature = surface_air_temperature
        self.dewpoint_temperature = dewpoint_temperature
        self.surface_downward_solar_radiation = surface_downward_solar_radiation
        self.surface_solar_radiation = surface_solar_radiation
        self.surface_thermal_radiation = surface_thermal_radiation
        self.surface_downward_thermal_radiation = surface_downward_thermal_radiation
        self.total_rain=total_rain
        self.ncfile = None
        self.update = update;

        if self.update == False :
            self.ncfile = Dataset(self.filename, 'w', format='NETCDF4')
        else:
            self.ncfile = Dataset(self.filename, 'r+', format='NETCDF4')

        self.ncfile.description = 'ECMWF Writer. Generated with pySpatialETL'

        if self.update == False:
            # dimensions
            self.ncfile.createDimension('time', None)
            self.ncfile.createDimension('lat', len(self.y_axis))
            self.ncfile.createDimension('lon', len(self.x_axis))

            # variables
            times = self.ncfile.createVariable('time', float64, ('time',))
            times.units= 'seconds since 2008-01-29 00:00:00'
            times.calendar= 'gregorian'
            times.standard_name= 'time'
            times.axis='T'
            times.conventions = "UTC time"

            latitudes = self.ncfile.createVariable('lat', float32, ('lat',))
            latitudes.units = "degree_north" ;
            latitudes.long_name = "latitude" ;
            latitudes.standard_name = "latitude" ;
            latitudes.valid_min = -90.0;
            latitudes.valid_max = 90.0 ;
            latitudes.axis = "Y" ;

            longitudes = self.ncfile.createVariable('lon', float32, ('lon',))
            longitudes.units = "degree_east" ;
            longitudes.long_name = "longitude" ;
            longitudes.standard_name = "longitude" ;
            longitudes.valid_min = -180.0 ;
            longitudes.valid_max = 180.0 ;
            longitudes.axis = "X" ;

            # data
            latitudes[:] = self.y_axis;
            longitudes[:] = self.x_axis;
            times[:] = date2num(self.t_axis,units = times.units, calendar = times.calendar);

    def close(self):
        self.ncfile.close()

    def write_variable_3D_mask(self):

        var = self.ncfile.createVariable('LSM', float32, ('time', 'lat', 'lon',), fill_value=9.96921e+36)
        var.long_name = "Land/sea mask";
        var.code = 172;
        var.table = 128;
        var[::] = 0

    def write_variable_surface_pressure(self):

        if self.update == False:
            var = self.ncfile.createVariable('SP', float32, ('time', 'lat', 'lon',), fill_value=9.96921e+36)
            var.long_name = "Surface pressure";
            var.code = 134;
            var.table = 128;
            var.units = "Pa";
        else:
            var = self.ncfile.variables['SP'];

        var[:] = self.surface_air_pressure*100;

    def write_variable_sea_surface_air_pressure(self):

        if self.update == False:
            var = self.ncfile.createVariable('MSL', float32, ('time', 'lat', 'lon',), fill_value=9.96921e+36)
            var.long_name = "Sea surface air pressure";
            var.code = 134;
            var.table = 128;
            var.units = "Pa";
        else:
            var = self.ncfile.variables['MSL'];

        var[:] = self.sea_surface_air_pressure*100;

    def write_variable_wind(self):

        if self.update == False:
            data_u = self.ncfile.createVariable('U10M', float32, ('time', 'lat', 'lon',), fill_value=9.96921e+36)
            data_u.long_name = "10 metre U wind component";
            data_u.code = 147;
            data_u.table = 128;
            data_u.units = "m s**-1";

            data_v = self.ncfile.createVariable('V10M', float32, ('time', 'lat', 'lon',), fill_value=9.96921e+36)
            data_v.long_name = "10 metre V wind component";
            data_v.code = 147;
            data_v.table = 128;
            data_v.units = "m s**-1";
        else :
            data_u = self.ncfile.variables['U10M']
            data_v = self.ncfile.variables['V10M']

        #data_u[:] = 270. - (math.degrees(self.wind_speed*math.sin(math.radians(self.wind_to_direction_angle)))) + 180.0 % 360.0
        #data_v[:] = 270. - (math.degrees(self.wind_speed*math.cos(math.radians(self.wind_to_direction_angle)) + 180.0 % 360.0

        data_u[:] = self.wind_speed * math.sin(math.radians((self.wind_from_direction_angle + 180.0) % 360.0))
        data_v[:] = self.wind_speed * math.cos(math.radians((self.wind_from_direction_angle + 180.0) % 360.0))

    def write_variable_surface_downward_sensible_heat_flux(self):

        if self.update == False:
            var = self.ncfile.createVariable('SSHF', float32, ('time', 'lat', 'lon',), fill_value=9.96921e+36)
            var.long_name = "Surface sensible heat flux";
            var.code = 146;
            var.table = 128;
            var.units = "W m**-2 s";
        else:
            var = self.ncfile.variables['SSHF']

        var[:] = self.surface_downward_sensible_heat_flux

    def write_variable_surface_downward_latent_heat_flux(self):

        var = self.ncfile.createVariable('SLHF', float32, ('time', 'lat', 'lon',), fill_value=9.96921e+36)
        var.long_name = "Surface latent heat flux";
        var.code = 147;
        var.table = 128;
        var.units = "W m**-2 s";
        var[:] = self.surface_downward_latent_heat_flux

    def write_variable_surface_air_temperature(self):

        var = self.ncfile.createVariable('T2M', float32, ('time', 'lat', 'lon',), fill_value=9.96921e+36)
        var.long_name = "2 metre temperature";
        var.code = 167;
        var.table = 128;
        var.units = "K";
        var[:] = self.surface_air_temperature

    def write_variable_dewpoint_temperature(self):

        var = self.ncfile.createVariable('D2M', float32, ('time', 'lat', 'lon',), fill_value=9.96921e+36)
        var.long_name = "2 metre dewpoint temperature";
        var.code = 168;
        var.table = 128;
        var.units = "K";
        # Value
        var[:] = self.dewpoint_temperature

    def write_variable_surface_downward_solar_radiation(self):

        var = self.ncfile.createVariable('SSRD', float32, ('time', 'lat', 'lon',), fill_value=9.96921e+36)
        var.long_name = "Surface solar radiation downwards";
        var.code = 169;
        var.table = 128;
        var.units = "W m**-2 s";
        var[:] = self.surface_downward_solar_radiation

    def write_variable_surface_downward_thermal_radiation(self):

        var = self.ncfile.createVariable('STRD', float32, ('time', 'lat', 'lon',), fill_value=9.96921e+36)
        var.long_name = "Surface thermal radiation downwards";
        var.code = 175;
        var.table = 128;
        var.units = "W m**-2 s";
        var[:] = self.surface_downward_thermal_radiation

    def write_variable_surface_solar_radiation(self):

        var = self.ncfile.createVariable('SSR', float32, ('time', 'lat', 'lon',), fill_value=9.96921e+36)
        var.long_name = "Surface solar radiation";
        var.code = 176;
        var.table = 128;
        var.units = "W m**-2 s";
        var[:] = self.surface_solar_radiation;

    def write_variable_surface_thermal_radiation(self):

        var = self.ncfile.createVariable('STR', float32, ('time', 'lat', 'lon',), fill_value=9.96921e+36)
        var.long_name = "Surface thermal radiation";
        var.code = 177;
        var.table = 128;
        var.units = "W m**-2 s";
        var[:] = self.surface_thermal_radiation,

    def write_variable_rain(self):

        var = self.ncfile.createVariable('TP', float32, ('time', 'lat', 'lon',), fill_value=9.96921e+36)
        var.long_name = "Total precipitation";
        var.code = 228;
        var.table = 128;
        var.units = "m";
        var[:] = self.total_rain