#! /usr/bin/env python2.7
# -*- coding: utf-8 -*-
#
# CoverageProcessing is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# any later version.
#
# CoverageProcessing is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# Author : Fabien Rétif - fabien.retif@zoho.com
#
from __future__ import division, print_function, absolute_import
from point.io.MultiPointWriter import MultiPointWriter
from point.TimeMultiPoint import TimeMultiPoint
from utils.VariableUnits import VariableUnits
from netCDF4 import Dataset
from netCDF4 import date2num
from numpy import float32,float64,int32

import numpy as np
import logging

class DefaultTimeMultiPointWriter(MultiPointWriter):

    def __init__(self, p,myFile):
        MultiPointWriter.__init__(self,p,myFile)

        if not isinstance(self.points, TimeMultiPoint):
            raise ValueError("This writer supports only a TimeMultiPoint object")

        self.ncfile = Dataset(self.filename, 'w', format='NETCDF4')
        self.ncfile.description = 'Generated with pyGeoSpatialETL'

        self.ncfile.data_source = str(self.points.data_source)
        self.ncfile.meta_data = str(self.points.meta_data)

        # Time dimension
        self.ncfile.createDimension('time', None)
        times = self.ncfile.createVariable('time', float64, ('time',))
        times.units= 'seconds since 1970-01-01 00:00:00'
        times.calendar= 'gregorian'
        times.standard_name= 'time'
        times.axis='T'
        times.conventions = "UTC time"

        times[:] = date2num(self.points.read_axis_t(), units = times.units, calendar = times.calendar)

        # Profil dimension
        self.ncfile.createDimension('point', self.points.get_nb_points())
        var = self.ncfile.createVariable('point', int32, ('point',))
        var.long_name = "Point number";
        var.standard_name = "point_number";
        var.axis = "X";
        # data
        var[:] = range(0,self.points.get_nb_points());

        var = self.ncfile.createVariable('latitude', float32, ('point',), fill_value=9.96921e+36)
        var.long_name = "latitude";
        var.standard_name = "latitude";
        var.valid_min = "-90.0";
        var.valid_max = "90.0";
        var.units = VariableUnits.CANONICAL_UNITS[var.standard_name];
        var[:] = self.points.read_axis_y()

        var = self.ncfile.createVariable('longitude', float32, ('point',), fill_value=9.96921e+36)
        var.long_name = "longitude";
        var.standard_name = "longitude";
        var.valid_min = "-180.0";
        var.valid_max = "180.0";
        var.units = VariableUnits.CANONICAL_UNITS[var.standard_name];
        var[:] = self.points.read_axis_x()


    def close(self):
        self.ncfile.close()

    def write_variable_wind_speed_10m(self):
        var = self.ncfile.createVariable('wind_speed_10m', float32, ('time','point',),
                                         fill_value=9.96921e+36)
        var.long_name = "Wind Speed 10m";
        var.standard_name = "wind_speed_10m";
        var.units = VariableUnits.CANONICAL_UNITS[var.standard_name];

        logging.info('[DefaultTimeMultiPointWriter] Writing variable \'' + var.long_name + '\'')

        time_index = 0
        for time in self.points.read_axis_t():
            logging.debug(
                '[DefaultTimeMultiPointWriter] Writing variable \'' + var.long_name + '\' at time \'' + str(time) + '\'')

            var[time_index:time_index + 1, :] = self.points.read_variable_wind_speed_10m_at_time(time)

            time_index += 1

    def write_variable_wind_from_direction_10m(self):
        var = self.ncfile.createVariable('wind_from_direction_10m', float32, ('time','point',),
                                         fill_value=9.96921e+36)
        var.long_name = "Wind From Direction 10m";
        var.standard_name = "wind_from_direction_10m";
        var.units = VariableUnits.CANONICAL_UNITS[var.standard_name];

        logging.info('[DefaultTimeMultiPointWriter] Writing variable \'' + var.long_name + '\'')

        time_index = 0
        for time in self.points.read_axis_t():
            logging.debug(
                '[DefaultTimeMultiPointWriter] Writing variable \'' + var.long_name + '\' at time \'' + str(time) + '\'')

            var[time_index:time_index + 1, :] = self.points.read_variable_wind_from_direction_10m_at_time(time)

            time_index += 1

    def write_variable_wind_to_direction_10m(self):
        var = self.ncfile.createVariable('wind_to_direction_10m', float32, ('time', 'point',),
                                         fill_value=9.96921e+36)
        var.long_name = "Wind To Direction 10m";
        var.standard_name = "wind_to_direction_10m";
        var.units = VariableUnits.CANONICAL_UNITS[var.standard_name];

        logging.info('[DefaultTimeMultiPointWriter] Writing variable \'' + var.long_name + '\'')

        time_index = 0
        for time in self.points.read_axis_t():
            logging.debug(
                '[DefaultTimeMultiPointWriter] Writing variable \'' + var.long_name + '\' at time \'' + str(
                    time) + '\'')

            var[time_index:time_index + 1, :] = self.points.read_variable_wind_to_direction_10m_at_time(time)

            time_index += 1

    def write_variable_rainfall_amount(self):
        var = self.ncfile.createVariable('rainfall_amount', float32, ('time','point',),
                                         fill_value=9.96921e+36)
        var.long_name = "Rainfall Amount";
        var.standard_name = "rainfall_amount";
        var.units = VariableUnits.CANONICAL_UNITS[var.standard_name];

        logging.info('[DefaultTimeMultiPointWriter] Writing variable \'' + var.long_name + '\'')

        time_index = 0
        for time in self.points.read_axis_t():
            logging.debug(
                '[DefaultTimeMultiPointWriter] Writing variable \'' + var.long_name + '\' at time \'' + str(time) + '\'')

            var[time_index:time_index + 1, :] = self.points.read_variable_rainfall_amount_at_time(time)

            time_index += 1

    def write_variable_surface_air_pressure(self):
        var = self.ncfile.createVariable('surface_air_pressure', float32, ('time','point'), fill_value=9.96921e+36)
        var.long_name = "Surface Air Pressure";
        var.standard_name = "surface_air_pressure";
        var.units = VariableUnits.CANONICAL_UNITS[var.standard_name];

        logging.info('[DefaultTimeMultiPointWriter] Writing variable \'' + var.long_name + '\'')

        time_index = 0
        for time in self.points.read_axis_t():
            logging.debug(
                '[DefaultTimeMultiPointWriter] Writing variable \'' + var.long_name + '\' at time \'' + str(time) + '\'')

            var[time_index:time_index + 1, :] = self.points.read_variable_surface_air_pressure_at_time(time)

            time_index += 1

    def write_variable_sea_surface_temperature(self):
        var = self.ncfile.createVariable('sea_surface_temperature', float32, ('time', 'point'), fill_value=9.96921e+36)
        var.long_name = "Sea Surface Temperature";
        var.standard_name = "sea_surface_temperature";
        var.units = VariableUnits.CANONICAL_UNITS[var.standard_name];

        logging.info('[DefaultTimeMultiPointWriter] Writing variable \'' + var.long_name + '\'')

        time_index = 0
        for time in self.points.read_axis_t():
            logging.debug(
                '[DefaultTimeMultiPointWriter] Writing variable \'' + var.long_name + '\' at time \'' + str(time) + '\'')

            var[time_index:time_index + 1, :] = self.points.read_variable_sea_surface_temperature_at_time(time)

            time_index += 1

    def write_variable_sea_surface_salinity(self):
        var = self.ncfile.createVariable('sea_surface_salinity', float32, ('time', 'point'), fill_value=9.96921e+36)
        var.long_name = "Sea Surface Salinity";
        var.standard_name = "sea_surface_salinity";
        var.units = VariableUnits.CANONICAL_UNITS[var.standard_name];

        logging.info('[DefaultTimeMultiPointWriter] Writing variable \'' + var.long_name + '\'')

        time_index = 0
        for time in self.points.read_axis_t():
            logging.debug(
                '[DefaultTimeMultiPointWriter] Writing variable \'' + var.long_name + '\' at time \'' + str(
                    time) + '\'')

            var[time_index:time_index + 1, :] = self.points.read_variable_sea_surface_salinity_at_time(time)

            time_index += 1


    def write_variable_barotropic_sea_water_speed(self):
        var = self.ncfile.createVariable('barotropic_sea_water_speed', float32, ('time','point',),
                                         fill_value=9.96921e+36)
        var.long_name = "Barotropic Sea Water Speed";
        var.standard_name = "barotropic_sea_water_speed";
        var.units = VariableUnits.CANONICAL_UNITS[var.standard_name];

        logging.info('[DefaultTimeMultiPointWriter] Writing variable \'' + var.long_name + '\'')

        time_index = 0
        for time in self.points.read_axis_t():
            logging.debug(
                '[DefaultTimeMultiPointWriter] Writing variable \'' + var.long_name + '\' at time \'' + str(time) + '\'')

            var[time_index:time_index + 1, :] = self.points.read_variable_barotropic_sea_water_speed_at_time(time)

            time_index += 1

    def write_variable_barotropic_sea_water_from_direction(self):
        var = self.ncfile.createVariable('barotropic_sea_water_from_direction', float32, ('time','point',),
                                         fill_value=9.96921e+36)
        var.long_name = "Barotropic Sea Water From Direction";
        var.standard_name = "barotropic_sea_water_from_direction";
        var.units = VariableUnits.CANONICAL_UNITS[var.standard_name];

        logging.info('[DefaultTimeMultiPointWriter] Writing variable \'' + var.long_name + '\'')

        time_index = 0
        for time in self.points.read_axis_t():
            logging.debug(
                '[DefaultTimeMultiPointWriter] Writing variable \'' + var.long_name + '\' at time \'' + str(time) + '\'')

            var[time_index:time_index + 1, :] = self.points.read_variable_barotropic_sea_water_from_direction_at_time(time)

            time_index += 1

    def write_variable_barotropic_sea_water_to_direction(self):
        var = self.ncfile.createVariable('barotropic_sea_water_to_direction', float32, ('time', 'point',),
                                         fill_value=9.96921e+36)
        var.long_name = "Barotropic Sea Water To Direction";
        var.standard_name = "barotropic_sea_water_to_direction";
        var.units = VariableUnits.CANONICAL_UNITS[var.standard_name];

        logging.info('[DefaultTimeMultiPointWriter] Writing variable \'' + var.long_name + '\'')

        time_index = 0
        for time in self.points.read_axis_t():
            logging.debug(
                '[DefaultTimeMultiPointWriter] Writing variable \'' + var.long_name + '\' at time \'' + str(
                    time) + '\'')

            var[time_index:time_index + 1, :] = self.points.read_variable_barotropic_sea_water_to_direction_at_time(time)

            time_index += 1