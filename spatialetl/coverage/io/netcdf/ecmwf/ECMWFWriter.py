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
from __future__ import division, print_function, absolute_import

from mpi4py import MPI
from netCDF4 import Dataset
from netCDF4 import date2num
from numpy import float32
from numpy import float64

from spatialetl.coverage.io.CoverageWriter import CoverageWriter
from spatialetl.utils.logger import logging


class ECMWFWriter (CoverageWriter):

    def __init__(self,cov,myFile):
        CoverageWriter.__init__(self,cov,myFile);

        self.ncfile = None
        format = 'NETCDF4_CLASSIC'

        if self.coverage.is_regular_grid()==False:
            raise IOError("This writer is specific to regular grid.")

        self.ncfile = Dataset(self.filename, 'w', parallel=True, comm=self.coverage.comm, info=MPI.Info(),
                              format=format)
        self.ncfile.description = 'ECMWF Writer. Generated with Coverage Processing tools'

        # dimensions
        self.ncfile.createDimension('time', self.coverage.get_t_size(type="target_global"))
        self.ncfile.createDimension('lat', self.coverage.get_y_size(type="target_global"))
        self.ncfile.createDimension('lon', self.coverage.get_x_size(type="target_global"))

        # variables
        times = self.ncfile.createVariable('time', float64, ('time',))
        times.units= 'seconds since 1970-01-01 00:00:00' 
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
        latitudes[self.coverage.map_mpi[self.coverage.rank]["dst_global_y"]] = self.coverage.read_axis_y()
        longitudes[self.coverage.map_mpi[self.coverage.rank]["dst_global_x"]] = self.coverage.read_axis_x()
        times[
            self.coverage.map_mpi[self.coverage.rank]["dst_global_t"]] = date2num(self.coverage.read_axis_t(),
                                                                                  units=times.units,
                                                                                  calendar=times.calendar)

    def close(self):
        self.ncfile.close()

    def write_variable_3D_land_binary_mask(self):
        var = self.ncfile.createVariable('LSM', float32, ('time', 'lat', 'lon',), fill_value=9.96921e+36)
        var.long_name = "Land/sea mask";
        var.code = 172;
        var.table = 128;

        time_index = 0
        for time in self.coverage.read_axis_t():
            logging.info('[ECMWFWriter] Writing variable \'3D mask\' at time \'' + str(time) + '\'')
            var[
            self.coverage.map_mpi[self.coverage.rank]["dst_global_t"].start + time_index:
            self.coverage.map_mpi[self.coverage.rank]["dst_global_t"].start + time_index + 1,
            self.coverage.map_mpi[self.coverage.rank]["dst_global_y"],
            self.coverage.map_mpi[self.coverage.rank]["dst_global_x"]
            ] = self.coverage.read_variable_2D_land_binary_mask_at_time(time)

            time_index += 1

    def write_variable_surface_air_pressure(self):
        var = self.ncfile.createVariable('SP', float32, ('time', 'lat', 'lon',), fill_value=9.96921e+36)
        var.long_name = "Surface pressure";
        var.code = 134;
        var.table = 128;
        var.units = "Pa";

        time_index = 0
        for time in self.coverage.read_axis_t():
            logging.info('[ECMWFWriter] Writing variable \'Surface pressure\' at time \'' + str(time) + '\'')
            data =  self.coverage.read_variable_surface_pressure_at_time(time)
            data[:] *= 100  # hPa to Pa
            var[
            self.coverage.map_mpi[self.coverage.rank]["dst_global_t"].start + time_index:
            self.coverage.map_mpi[self.coverage.rank]["dst_global_t"].start + time_index + 1,
            self.coverage.map_mpi[self.coverage.rank]["dst_global_y"],
            self.coverage.map_mpi[self.coverage.rank]["dst_global_x"]
            ] = data
            time_index += 1

    def write_variable_sea_surface_air_pressure(self):
        var = self.ncfile.createVariable('MSL', float32, ('time', 'lat', 'lon',), fill_value=9.96921e+36)
        var.long_name = "Surface pressure";
        var.code = 134;
        var.table = 128;
        var.units = "Pa";

        time_index = 0
        for time in self.coverage.read_axis_t():
            logging.info('[ECMWFWriter] Writing variable \'Surface pressure\' at time \'' + str(time) + '\'')
            data = self.coverage.read_variable_sea_surface_air_pressure_at_time(time)
            data[:] *= 100 # hPa to Pa
            var[
            self.coverage.map_mpi[self.coverage.rank]["dst_global_t"].start + time_index:
            self.coverage.map_mpi[self.coverage.rank]["dst_global_t"].start + time_index + 1,
            self.coverage.map_mpi[self.coverage.rank]["dst_global_y"],
            self.coverage.map_mpi[self.coverage.rank]["dst_global_x"]
            ] = data
            time_index += 1

    def write_variable_surface_downward_sensible_heat_flux(self):
        var = self.ncfile.createVariable('SSHF', float32, ('time', 'lat', 'lon',), fill_value=9.96921e+36)
        var.long_name = "Surface sensible heat flux";
        var.code = 146;
        var.table = 128;
        var.units = "W m**-2 s";

        time_index = 0
        for time in self.coverage.read_axis_t():
            logging.info('[ECMWFWriter] Writing variable \'Surface sensible heat flux\' at time \'' + str(time) + '\'')

            var[
            self.coverage.map_mpi[self.coverage.rank]["dst_global_t"].start + time_index:
            self.coverage.map_mpi[self.coverage.rank]["dst_global_t"].start + time_index + 1,
            self.coverage.map_mpi[self.coverage.rank]["dst_global_y"],
            self.coverage.map_mpi[self.coverage.rank]["dst_global_x"]
            ] = self.coverage.read_variable_surface_downward_sensible_heat_flux_at_time(time)
            time_index += 1

    def write_variable_surface_downward_latent_heat_flux(self):
        var = self.ncfile.createVariable('SLHF', float32, ('time', 'lat', 'lon',), fill_value=9.96921e+36)
        var.long_name = "Surface latent heat flux";
        var.code = 147;
        var.table = 128;
        var.units = "W m**-2 s";

        time_index = 0
        for time in self.coverage.read_axis_t():
            logging.info('[ECMWFWriter] Writing variable \'Surface latent heat flux\' at time \'' + str(time) + '\'')
            var[
            self.coverage.map_mpi[self.coverage.rank]["dst_global_t"].start + time_index:
            self.coverage.map_mpi[self.coverage.rank]["dst_global_t"].start + time_index + 1,
            self.coverage.map_mpi[self.coverage.rank]["dst_global_y"],
            self.coverage.map_mpi[self.coverage.rank]["dst_global_x"]
            ] = self.coverage.read_variable_surface_downward_latent_heat_flux_at_time(time)
            time_index += 1

    def write_variable_wind_10m(self):
        ucomp = self.ncfile.createVariable('U10M', float32, ('time','lat', 'lon',),fill_value=9.96921e+36)
        ucomp.long_name = "10 metre U wind component" ;
        ucomp.code = 147;
        ucomp.table = 128;
        ucomp.units = "m s**-1" ;

        vcomp = self.ncfile.createVariable('V10M', float32, ('time', 'lat', 'lon',),fill_value=9.96921e+36)
        vcomp.long_name = "10 metre V wind component" ;
        vcomp.code = 147;
        vcomp.table = 128;
        vcomp.units = "m s**-1" ;


        time_index=0
        for time in self.coverage.read_axis_t():

            logging.info('[ECMWFWriter] Writing variable \'Wind\' at time \''+str(time)+'\'')

            data_u, data_v = self.coverage.read_variable_wind_10m_at_time(time)

            ucomp[
            self.coverage.map_mpi[self.coverage.rank]["dst_global_t"].start + time_index:
            self.coverage.map_mpi[self.coverage.rank]["dst_global_t"].start + time_index + 1,
            self.coverage.map_mpi[self.coverage.rank]["dst_global_y"],
            self.coverage.map_mpi[self.coverage.rank]["dst_global_x"]
            ] = data_u

            vcomp[
            self.coverage.map_mpi[self.coverage.rank]["dst_global_t"].start + time_index:
            self.coverage.map_mpi[self.coverage.rank]["dst_global_t"].start + time_index + 1,
            self.coverage.map_mpi[self.coverage.rank]["dst_global_y"],
            self.coverage.map_mpi[self.coverage.rank]["dst_global_x"]
            ] = data_v

            time_index += 1

    def write_variable_surface_air_temperature(self):
        var = self.ncfile.createVariable('T2M', float32, ('time', 'lat', 'lon',), fill_value=9.96921e+36)
        var.long_name = "2 metre temperature";
        var.code = 167;
        var.table = 128;
        var.units = "K";

        time_index = 0
        for time in self.coverage.read_axis_t():
            logging.info('[ECMWFWriter] Writing variable \'Surface air temperature\' at time \'' + str(time) + '\'')
            data = self.coverage.read_variable_surface_air_temperature_at_time(time)
            data += 273.15 # Celsius to Kelvin
            var[
            self.coverage.map_mpi[self.coverage.rank]["dst_global_t"].start + time_index:
            self.coverage.map_mpi[self.coverage.rank]["dst_global_t"].start + time_index + 1,
            self.coverage.map_mpi[self.coverage.rank]["dst_global_y"],
            self.coverage.map_mpi[self.coverage.rank]["dst_global_x"]
            ] = data
            time_index += 1

    def write_variable_dew_point_temperature(self):
        var = self.ncfile.createVariable('D2M', float32, ('time', 'lat', 'lon',), fill_value=9.96921e+36)
        var.long_name = "2 metre dewpoint temperature";
        var.code = 168;
        var.table = 128;
        var.units = "K";

        time_index = 0
        for time in self.coverage.read_axis_t():
            logging.info('[ECMWFWriter] Writing variable \'Dewpoint temperature\' at time \'' + str(time) + '\'')
            data = self.coverage.read_variable_dew_point_temperature_at_time(time)
            data += 273.15  # Celsius to Kelvin
            var[
            self.coverage.map_mpi[self.coverage.rank]["dst_global_t"].start + time_index:
            self.coverage.map_mpi[self.coverage.rank]["dst_global_t"].start + time_index + 1,
            self.coverage.map_mpi[self.coverage.rank]["dst_global_y"],
            self.coverage.map_mpi[self.coverage.rank]["dst_global_x"]
            ] = data
            time_index += 1

    def write_variable_surface_downward_solar_radiation(self):
        var = self.ncfile.createVariable('SSRD', float32, ('time', 'lat', 'lon',), fill_value=9.96921e+36)
        var.long_name = "Surface solar radiation downwards";
        var.code = 169;
        var.table = 128;
        var.units = "W m**-2 s";

        time_index = 0
        for time in self.coverage.read_axis_t():
            logging.info('[ECMWFWriter] Writing variable \'Surface solar radiation downwards\' at time \'' + str(time) + '\'')
            var[
            self.coverage.map_mpi[self.coverage.rank]["dst_global_t"].start + time_index:
            self.coverage.map_mpi[self.coverage.rank]["dst_global_t"].start + time_index + 1,
            self.coverage.map_mpi[self.coverage.rank]["dst_global_y"],
            self.coverage.map_mpi[self.coverage.rank]["dst_global_x"]
            ] = self.coverage.read_variable_surface_downward_solar_radiation_at_time(time)
            time_index += 1

    def write_variable_surface_downward_thermal_radiation(self):
        var = self.ncfile.createVariable('STRD', float32, ('time', 'lat', 'lon',), fill_value=9.96921e+36)
        var.long_name = "Surface thermal radiation downwards";
        var.code = 175;
        var.table = 128;
        var.units = "W m**-2 s";

        time_index = 0
        for time in self.coverage.read_axis_t():
            logging.info('[ECMWFWriter] Writing variable \'Surface thermal radiation downwards\' at time \'' + str(time) + '\'')
            var[
            self.coverage.map_mpi[self.coverage.rank]["dst_global_t"].start + time_index:
            self.coverage.map_mpi[self.coverage.rank]["dst_global_t"].start + time_index + 1,
            self.coverage.map_mpi[self.coverage.rank]["dst_global_y"],
            self.coverage.map_mpi[self.coverage.rank]["dst_global_x"]
            ] = self.coverage.read_variable_surface_downward_thermal_radiation_at_time(time)
            time_index += 1

    def write_variable_surface_solar_radiation(self):
        var = self.ncfile.createVariable('SSR', float32, ('time', 'lat', 'lon',), fill_value=9.96921e+36)
        var.long_name = "Surface solar radiation";
        var.code = 176;
        var.table = 128;
        var.units = "W m**-2 s";

        time_index = 0
        for time in self.coverage.read_axis_t():
            logging.info('[ECMWFWriter] Writing variable \'Surface solar radiation\' at time \'' + str(time) + '\'')
            var[
            self.coverage.map_mpi[self.coverage.rank]["dst_global_t"].start + time_index:
            self.coverage.map_mpi[self.coverage.rank]["dst_global_t"].start + time_index + 1,
            self.coverage.map_mpi[self.coverage.rank]["dst_global_y"],
            self.coverage.map_mpi[self.coverage.rank]["dst_global_x"]
            ] = self.coverage.read_variable_surface_solar_radiation_at_time(time)

            time_index += 1

    def write_variable_surface_thermal_radiation(self):
        var = self.ncfile.createVariable('STR', float32, ('time', 'lat', 'lon',), fill_value=9.96921e+36)
        var.long_name = "Surface thermal radiation";
        var.code = 177;
        var.table = 128;
        var.units = "W m**-2 s";

        time_index = 0
        for time in self.coverage.read_axis_t():
            logging.info('[ECMWFWriter] Writing variable \'Surface thermal radiation\' at time \'' + str(time) + '\'')
            var[
            self.coverage.map_mpi[self.coverage.rank]["dst_global_t"].start + time_index:
            self.coverage.map_mpi[self.coverage.rank]["dst_global_t"].start + time_index + 1,
            self.coverage.map_mpi[self.coverage.rank]["dst_global_y"],
            self.coverage.map_mpi[self.coverage.rank]["dst_global_x"]
            ] = self.coverage.read_variable_surface_thermal_radiation_at_time(time)
            time_index += 1

    def write_variable_rainfall_amount(self):
        var = self.ncfile.createVariable('TP', float32, ('time', 'lat', 'lon',), fill_value=9.96921e+36)
        var.long_name = "Total precipitation";
        var.code = 228;
        var.table = 128;
        var.units = "m";

        time_index = 0
        for time in self.coverage.read_axis_t():
            logging.info('[ECMWFWriter] Writing variable \'Rain\' at time \'' + str(time) + '\'')
            var[
            self.coverage.map_mpi[self.coverage.rank]["dst_global_t"].start + time_index:
            self.coverage.map_mpi[self.coverage.rank]["dst_global_t"].start + time_index + 1,
            self.coverage.map_mpi[self.coverage.rank]["dst_global_y"],
            self.coverage.map_mpi[self.coverage.rank]["dst_global_x"]
            ] = self.coverage.read_variable_rainfall_amount_at_time(time)
            time_index += 1