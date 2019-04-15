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
from utils.VariableDefinition import VariableDefinition
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

        # Geo-points dimension
        self.ncfile.createDimension(VariableDefinition.VARIABLE_NAME['point'], self.points.get_nb_points())
        var = self.ncfile.createVariable(VariableDefinition.VARIABLE_NAME['point'], int32,
                                         (VariableDefinition.VARIABLE_NAME['point'],))
        var.long_name = VariableDefinition.LONG_NAME['point']
        var.standard_name = VariableDefinition.STANDARD_NAME['point']
        var.axis = "X";
        var.units = VariableDefinition.CANONICAL_UNITS['point'];
        var[:] = range(0, self.points.get_nb_points());

        var = self.ncfile.createVariable(VariableDefinition.VARIABLE_NAME['latitude'], float32,
                                         (VariableDefinition.VARIABLE_NAME['point'],), fill_value=9.96921e+36)
        var.long_name = VariableDefinition.LONG_NAME['latitude']
        var.standard_name = VariableDefinition.STANDARD_NAME['latitude']
        var.valid_min = "-90.0";
        var.valid_max = "90.0";
        var.units = VariableDefinition.CANONICAL_UNITS['latitude']
        var[:] = self.points.read_axis_y()

        var = self.ncfile.createVariable(VariableDefinition.VARIABLE_NAME['longitude'], float32,
                                         (VariableDefinition.VARIABLE_NAME['point'],), fill_value=9.96921e+36)
        var.long_name = VariableDefinition.LONG_NAME['longitude']
        var.standard_name = VariableDefinition.STANDARD_NAME['longitude']
        var.valid_min = "-180.0";
        var.valid_max = "180.0";
        var.units = VariableDefinition.CANONICAL_UNITS['longitude']
        var[:] = self.points.read_axis_x()

        # Time dimension
        self.ncfile.createDimension(VariableDefinition.VARIABLE_NAME['time'], None)
        times = self.ncfile.createVariable(VariableDefinition.VARIABLE_NAME['time'], float64, (VariableDefinition.VARIABLE_NAME['time'],))
        times.units = 'seconds since 1970-01-01 00:00:00'
        times.calendar = 'gregorian'
        times.standard_name = 'time'
        times.axis = 'T'
        times.conventions = "UTC time"

        times[:] = date2num(self.points.read_axis_t(), units=times.units, calendar=times.calendar)


    def close(self):
        self.ncfile.close()

    def write_variable_wind_speed_10m(self):
        var = self.ncfile.createVariable(VariableDefinition.VARIABLE_NAME['wind_speed_10m'], float32, (VariableDefinition.VARIABLE_NAME['time'],VariableDefinition.VARIABLE_NAME['point'],),
                                         fill_value=9.96921e+36)
        var.long_name = VariableDefinition.LONG_NAME['wind_speed_10m']
        var.standard_name = VariableDefinition.STANDARD_NAME['wind_speed_10m']
        var.units = VariableDefinition.CANONICAL_UNITS['wind_speed_10m']

        logging.info('[DefaultTimeMultiPointWriter] Writing variable \'' + str(VariableDefinition.LONG_NAME['wind_speed_10m']) + '\'')

        time_index = 0
        for time in self.points.read_axis_t():
            logging.debug(
                '[DefaultTimeMultiPointWriter] Writing variable \'' + str(VariableDefinition.LONG_NAME['wind_speed_10m']) + '\' at time \'' + str(time) + '\'')

            var[time_index:time_index + 1, :] = self.points.read_variable_wind_speed_10m_at_time(time)

            time_index += 1

    def write_variable_wind_from_direction_10m(self):
        var = self.ncfile.createVariable(VariableDefinition.VARIABLE_NAME['wind_from_direction_10m'], float32, (VariableDefinition.VARIABLE_NAME['time'],VariableDefinition.VARIABLE_NAME['point'],),
                                         fill_value=9.96921e+36)
        var.long_name = VariableDefinition.LONG_NAME['wind_from_direction_10m']
        var.standard_name = VariableDefinition.STANDARD_NAME['wind_from_direction_10m']
        var.units = VariableDefinition.CANONICAL_UNITS['wind_from_direction_10m']

        logging.info('[DefaultTimeMultiPointWriter] Writing variable \'' + str(VariableDefinition.LONG_NAME['wind_from_direction_10m']) + '\'')

        time_index = 0
        for time in self.points.read_axis_t():
            logging.debug(
                '[DefaultTimeMultiPointWriter] Writing variable \'' + str(VariableDefinition.LONG_NAME['wind_from_direction_10m']) + '\' at time \'' + str(time) + '\'')

            var[time_index:time_index + 1, :] = self.points.read_variable_wind_from_direction_10m_at_time(time)

            time_index += 1

    def write_variable_wind_to_direction_10m(self):
        var = self.ncfile.createVariable(VariableDefinition.VARIABLE_NAME['wind_to_direction_10m'], float32, (VariableDefinition.VARIABLE_NAME['time'],VariableDefinition.VARIABLE_NAME['point'],),
                                         fill_value=9.96921e+36)
        var.long_name = VariableDefinition.LONG_NAME['wind_to_direction_10m']
        var.standard_name = VariableDefinition.STANDARD_NAME['wind_to_direction_10m']
        var.units = VariableDefinition.CANONICAL_UNITS['wind_to_direction_10m']

        logging.info('[DefaultTimeMultiPointWriter] Writing variable \'' + str(VariableDefinition.LONG_NAME['wind_to_direction_10m']) + '\'')

        time_index = 0
        for time in self.points.read_axis_t():
            logging.debug(
                '[DefaultTimeMultiPointWriter] Writing variable \'' + str(VariableDefinition.LONG_NAME['wind_to_direction_10m'])+'\' at time \'' + str(
                    time) + '\'')

            var[time_index:time_index + 1, :] = self.points.read_variable_wind_to_direction_10m_at_time(time)

            time_index += 1

    def write_variable_rainfall_amount(self):
        var = self.ncfile.createVariable(VariableDefinition.VARIABLE_NAME['rainfall_amount'], float32, (VariableDefinition.VARIABLE_NAME['time'],VariableDefinition.VARIABLE_NAME['point'],),
                                         fill_value=9.96921e+36)
        var.long_name = VariableDefinition.LONG_NAME['rainfall_amount']
        var.standard_name = VariableDefinition.STANDARD_NAME['rainfall_amount']
        var.units = VariableDefinition.CANONICAL_UNITS['rainfall_amount']

        logging.info('[DefaultTimeMultiPointWriter] Writing variable \'' + str(VariableDefinition.LONG_NAME['rainfall_amount']) + '\'')

        time_index = 0
        for time in self.points.read_axis_t():
            logging.debug(
                '[DefaultTimeMultiPointWriter] Writing variable \'' + str(VariableDefinition.LONG_NAME['rainfall_amount']) + '\' at time \'' + str(time) + '\'')

            var[time_index:time_index + 1, :] = self.points.read_variable_rainfall_amount_at_time(time)

            time_index += 1

    def write_variable_surface_air_pressure(self):
        var = self.ncfile.createVariable(VariableDefinition.VARIABLE_NAME['surface_air_pressure'], float32, (VariableDefinition.VARIABLE_NAME['time'],VariableDefinition.VARIABLE_NAME['point']), fill_value=9.96921e+36)
        var.long_name = VariableDefinition.LONG_NAME['surface_air_pressure']
        var.standard_name = VariableDefinition.STANDARD_NAME['surface_air_pressure']
        var.units = VariableDefinition.CANONICAL_UNITS['surface_air_pressure']

        logging.info('[DefaultTimeMultiPointWriter] Writing variable \'' + str(VariableDefinition.LONG_NAME['surface_air_pressure'])+ '\'')

        time_index = 0
        for time in self.points.read_axis_t():
            logging.debug(
                '[DefaultTimeMultiPointWriter] Writing variable \'' + str(VariableDefinition.LONG_NAME['surface_air_pressure']) + '\' at time \'' + str(time) + '\'')

            var[time_index:time_index + 1, :] = self.points.read_variable_surface_air_pressure_at_time(time)

            time_index += 1

    def write_variable_sea_surface_temperature(self):
        var = self.ncfile.createVariable(VariableDefinition.VARIABLE_NAME['sea_surface_temperature'], float32, (VariableDefinition.VARIABLE_NAME['time'],VariableDefinition.VARIABLE_NAME['point']), fill_value=9.96921e+36)
        var.long_name = VariableDefinition.LONG_NAME['sea_surface_temperature'];
        var.standard_name =VariableDefinition.STANDARD_NAME['sea_surface_temperature']
        var.units = VariableDefinition.CANONICAL_UNITS['sea_surface_temperature']

        logging.info('[DefaultTimeMultiPointWriter] Writing variable \'' + str(VariableDefinition.LONG_NAME['sea_surface_temperature']) + '\'')

        time_index = 0
        for time in self.points.read_axis_t():
            logging.debug(
                '[DefaultTimeMultiPointWriter] Writing variable \'' + str(VariableDefinition.LONG_NAME['sea_surface_temperature']) + '\' at time \'' + str(time) + '\'')

            var[time_index:time_index + 1, :] = self.points.read_variable_sea_surface_temperature_at_time(time)

            time_index += 1

    def write_variable_sea_surface_salinity(self):
        var = self.ncfile.createVariable(VariableDefinition.VARIABLE_NAME['sea_surface_salinity'], float32, (VariableDefinition.VARIABLE_NAME['time'],VariableDefinition.VARIABLE_NAME['point']), fill_value=9.96921e+36)
        var.long_name = VariableDefinition.LONG_NAME['sea_surface_salinity']
        var.standard_name = VariableDefinition.STANDARD_NAME['sea_surface_salinity']
        var.units = VariableDefinition.CANONICAL_UNITS['sea_surface_salinity']

        logging.info('[DefaultTimeMultiPointWriter] Writing variable \'' + str(VariableDefinition.LONG_NAME['sea_surface_salinity']) + '\'')

        time_index = 0
        for time in self.points.read_axis_t():
            logging.debug(
                '[DefaultTimeMultiPointWriter] Writing variable \'' + str(VariableDefinition.LONG_NAME['sea_surface_salinity']) + '\' at time \'' + str(
                    time) + '\'')

            var[time_index:time_index + 1, :] = self.points.read_variable_sea_surface_salinity_at_time(time)

            time_index += 1


    def write_variable_barotropic_sea_water_speed(self):
        var = self.ncfile.createVariable(VariableDefinition.VARIABLE_NAME['barotropic_sea_water_speed'], float32, (VariableDefinition.VARIABLE_NAME['time'],VariableDefinition.VARIABLE_NAME['point']),
                                         fill_value=9.96921e+36)
        var.long_name = VariableDefinition.LONG_NAME['barotropic_sea_water_speed']
        var.standard_name = VariableDefinition.STANDARD_NAME['barotropic_sea_water_speed']
        var.units = VariableDefinition.CANONICAL_UNITS['barotropic_sea_water_speed']

        logging.info('[DefaultTimeMultiPointWriter] Writing variable \'' + str(VariableDefinition.LONG_NAME['barotropic_sea_water_speed']) + '\'')

        time_index = 0
        for time in self.points.read_axis_t():
            logging.debug(
                '[DefaultTimeMultiPointWriter] Writing variable \'' + str(VariableDefinition.LONG_NAME['barotropic_sea_water_speed']) + '\' at time \'' + str(time) + '\'')

            var[time_index:time_index + 1, :] = self.points.read_variable_barotropic_sea_water_speed_at_time(time)

            time_index += 1

    def write_variable_barotropic_sea_water_from_direction(self):
        var = self.ncfile.createVariable(VariableDefinition.VARIABLE_NAME['barotropic_sea_water_from_direction'], float32, (VariableDefinition.VARIABLE_NAME['time'],VariableDefinition.VARIABLE_NAME['point']),
                                         fill_value=9.96921e+36)
        var.long_name = VariableDefinition.LONG_NAME['barotropic_sea_water_from_direction']
        var.standard_name = VariableDefinition.STANDARD_NAME['barotropic_sea_water_from_direction']
        var.units = VariableDefinition.CANONICAL_UNITS['barotropic_sea_water_from_direction']

        logging.info('[DefaultTimeMultiPointWriter] Writing variable \'' + str(VariableDefinition.LONG_NAME['barotropic_sea_water_from_direction']) + '\'')

        time_index = 0
        for time in self.points.read_axis_t():
            logging.debug(
                '[DefaultTimeMultiPointWriter] Writing variable \'' + str(VariableDefinition.LONG_NAME['barotropic_sea_water_from_direction']) + '\' at time \'' + str(time) + '\'')

            var[time_index:time_index + 1, :] = self.points.read_variable_barotropic_sea_water_from_direction_at_time(time)

            time_index += 1

    def write_variable_barotropic_sea_water_to_direction(self):
        var = self.ncfile.createVariable(VariableDefinition.VARIABLE_NAME['barotropic_sea_water_to_direction'], float32, (VariableDefinition.VARIABLE_NAME['time'],VariableDefinition.VARIABLE_NAME['point']),
                                         fill_value=9.96921e+36)
        var.long_name = VariableDefinition.LONG_NAME['barotropic_sea_water_to_direction']
        var.standard_name = VariableDefinition.STANDARD_NAME['barotropic_sea_water_to_direction']
        var.units = VariableDefinition.CANONICAL_UNITS['barotropic_sea_water_to_direction']

        logging.info('[DefaultTimeMultiPointWriter] Writing variable \'' + str(VariableDefinition.LONG_NAME['barotropic_sea_water_to_direction']) + '\'')

        time_index = 0
        for time in self.points.read_axis_t():
            logging.debug(
                '[DefaultTimeMultiPointWriter] Writing variable \'' + str(VariableDefinition.LONG_NAME['barotropic_sea_water_to_direction']) + '\' at time \'' + str(
                    time) + '\'')

            var[time_index:time_index + 1, :] = self.points.read_variable_barotropic_sea_water_to_direction_at_time(time)

            time_index += 1

    def write_variable_sea_water_pressure_at_sea_water_surface(self):
        var = self.ncfile.createVariable(VariableDefinition.VARIABLE_NAME['sea_water_pressure_at_sea_water_surface'], float32, (
        VariableDefinition.VARIABLE_NAME['time'], VariableDefinition.VARIABLE_NAME['point']), fill_value=9.96921e+36)
        var.long_name = VariableDefinition.LONG_NAME['sea_water_pressure_at_sea_water_surface']
        var.standard_name = VariableDefinition.STANDARD_NAME['sea_water_pressure_at_sea_water_surface']
        var.units = VariableDefinition.CANONICAL_UNITS['sea_water_pressure_at_sea_water_surface']

        logging.info('[DefaultTimeMultiPointWriter] Writing variable \'' + str(
            VariableDefinition.LONG_NAME['sea_water_pressure_at_sea_water_surface']) + '\'')

        time_index = 0
        for time in self.points.read_axis_t():
            logging.debug(
                '[DefaultTimeMultiPointWriter] Writing variable \'' + str(
                    VariableDefinition.LONG_NAME['sea_water_pressure_at_sea_water_surface']) + '\' at time \'' + str(
                    time) + '\'')

            var[time_index:time_index + 1, :] = self.points.read_variable_sea_water_pressure_at_sea_water_surface_at_time(time)

            time_index += 1

    def write_variable_water_volume_transport_into_sea_water_from_rivers(self):
        var = self.ncfile.createVariable(VariableDefinition.VARIABLE_NAME['water_volume_transport_into_sea_water_from_rivers'], float32, (
        VariableDefinition.VARIABLE_NAME['time'], VariableDefinition.VARIABLE_NAME['point']),
                                         fill_value=9.96921e+36)
        var.long_name = VariableDefinition.LONG_NAME['water_volume_transport_into_sea_water_from_rivers']
        var.standard_name = VariableDefinition.STANDARD_NAME['water_volume_transport_into_sea_water_from_rivers']
        var.units = VariableDefinition.CANONICAL_UNITS['water_volume_transport_into_sea_water_from_rivers']

        logging.info('[DefaultTimeMultiPointWriter] Writing variable \'' + str(
            VariableDefinition.LONG_NAME['water_volume_transport_into_sea_water_from_rivers']) + '\'')

        time_index = 0
        for time in self.points.read_axis_t():
            logging.debug(
                '[DefaultTimeMultiPointWriter] Writing variable \'' + str(
                    VariableDefinition.LONG_NAME['water_volume_transport_into_sea_water_from_rivers']) + '\' at time \'' + str(time) + '\'')

            var[time_index:time_index + 1, :] = self.points.read_variable_water_volume_transport_into_sea_water_from_rivers_at_time(time)

            time_index += 1