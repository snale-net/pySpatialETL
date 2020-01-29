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
from spatialetl.point.io.MultiPointWriter import MultiPointWriter
from spatialetl.point.TimeMultiPoint import TimeMultiPoint
from spatialetl.utils.VariableDefinition import VariableDefinition
from netCDF4 import Dataset
from netCDF4 import date2num
from numpy import float32
from numpy import float64
import numpy as np
import logging
import pandas


class DefaultTimePointWriter(MultiPointWriter):

    def __init__(self, s,index_point,myFile):
        MultiPointWriter.__init__(self, s, myFile);

        if not isinstance(self.points, TimeMultiPoint):
            raise ValueError("This writer supports only TimeMultiPoint object")

        self.index_x = index_point

        index = pandas.DatetimeIndex(self.points.read_axis_t())
        self.data = pandas.DataFrame(index=index)

    def close(self):
        self.data.to_csv(self.filename, sep='\t', columns=list(self.data), header=False, encoding='utf-8', na_rep="NaN")

        file = open(self.filename, "r+")
        old = file.read()  # read everything in the file
        file.seek(0)  # rewind

        file.write("############################################################ \n\
# Station : " + str(self.points.read_variable_point_names()[self.index_x]) + " \n\
# Coordinate Reference System : WGS84 \n\
# Longitude : " + str(self.points.read_axis_x()[self.index_x]) + " \n\
# Latitude : " + str(self.points.read_axis_y()[self.index_x]) + " \n\
# Data source : " + str(self.points.data_source) + " \n\
# Meta Data : " + str(self.points.meta_data) + " \n\
# Time zone : UTC \n\
# Separator: Tabulation \\t \n\
# Column 1: year-month-day hour:minute:second UTC \n")

        column = 2
        for key in list(self.data):
            file.write("# Column " + str(column) + ": " + str(VariableDefinition.STANDARD_NAME[key]) + " (" + str(
                VariableDefinition.CANONICAL_UNITS[key]) + ") - FillValue: NaN \n")
            column = column + 1

        file.write("############################################################\n")

        file.write(old)  # write the new line before
        file.close()

    def write_variable_sea_surface_height_above_mean_sea_level(self):
        logging.info('[DefaultTimePointWriter] Writing variable \''+str(VariableDefinition.LONG_NAME['sea_surface_height_above_mean_sea_level'])+'\'')

        data = np.zeros([self.points.get_t_size()])
        data[:] = np.nan
        time_index = 0
        for time in self.points.read_axis_t():
            logging.debug(
                '[DefaultTimePointWriter] Writing variable \'' + str(
                    VariableDefinition.LONG_NAME['sea_surface_height_above_mean_sea_level']) + '\' at time \'' + str(
                    time) + '\'')

            if self.points.is_raw_times:
                # Inutile de rechercher la date si l'axe du temps est le même que la donnée brute
                data[time_index] = self.points.read_variable_sea_surface_height_above_mean_sea_level_at_time(time_index)[self.index_x]
            else:
                data[time_index] = self.points.read_variable_sea_surface_height_above_mean_sea_level_at_time(time)[self.index_x]

            time_index += 1

        self.data['sea_surface_height_above_mean_sea_level'] = data

    def write_variable_sea_surface_height_above_geoid(self):
        logging.info('[DefaultTimePointWriter] Writing variable \'' + str(
            VariableDefinition.LONG_NAME['sea_surface_height_above_geoid']) + '\'')

        data = np.zeros([self.points.get_t_size()])
        data[:] = np.nan
        time_index = 0
        for time in self.points.read_axis_t():

            logging.debug(
                '[DefaultTimePointWriter] Writing variable \'' + str(
                    VariableDefinition.LONG_NAME['sea_surface_height_above_geoid']) + '\' at time \'' + str(
                    time) + '\'')

            if self.points.is_raw_times:
                # Inutile de rechercher la date si l'axe du temps est le même que la donnée brute
                data[time_index] = self.points.read_variable_sea_surface_height_at_time(time_index)[self.index_x]
            else:
                data[time_index] = self.points.read_variable_sea_surface_height_at_time(time)[self.index_x]

            time_index += 1

        self.data['sea_surface_height_above_geoid'] = data

    def write_variable_sea_surface_temperature(self):
        logging.info('[DefaultTimePointWriter] Writing variable \''+str(VariableDefinition.LONG_NAME['sea_surface_temperature'])+'\'')
        data = np.zeros([self.points.get_t_size()])
        data[:] = np.nan
        time_index = 0
        for time in self.points.read_axis_t():

            logging.debug(
                '[DefaultTimePointWriter] Writing variable \'' + str(
                    VariableDefinition.LONG_NAME['sea_surface_temperature']) + '\' at time \'' + str(
                    time) + '\'')

            if self.points.is_raw_times:
                # Inutile de rechercher la date si l'axe du temps est le même que la donnée brute
                data[time_index] = self.points.read_variable_sea_surface_temperature_at_time(time_index)[self.index_x]
            else:
                data[time_index] = self.points.read_variable_sea_surface_temperature_at_time(time)[self.index_x]

            time_index += 1

        self.data['sea_surface_temperature'] = data

    def write_variable_sea_water_temperature_at_ground_level(self):
        logging.info('[DefaultTimePointWriter] Writing variable \''+str(VariableDefinition.LONG_NAME['sea_water_temperature_at_ground_level'])+'\'')
        data = np.zeros([self.points.get_t_size()])
        data[:] = np.nan
        time_index = 0
        for time in self.points.read_axis_t():

            logging.debug(
                '[DefaultTimePointWriter] Writing variable \'' + str(
                    VariableDefinition.LONG_NAME['sea_water_temperature_at_ground_level']) + '\' at time \'' + str(
                    time) + '\'')

            if self.points.is_raw_times:
                # Inutile de rechercher la date si l'axe du temps est le même que la donnée brute
                data[time_index] = self.points.read_variable_sea_water_temperature_at_ground_level_at_time(time_index)[
                    self.index_x]
            else:
                data[time_index] = self.points.read_variable_sea_water_temperature_at_ground_level_at_time(time)[self.index_x]

            time_index += 1

        self.data['sea_water_temperature_at_ground_level'] = data

    def write_variable_longitude(self):
        logging.info('[DefaultTimePointWriter] Writing variable \''+str(VariableDefinition.LONG_NAME['longitude'])+'\'')
        data = np.zeros([self.points.get_t_size()])
        data[:] = np.nan
        time_index = 0
        for time in self.points.read_axis_t():

            logging.debug(
                '[DefaultTimePointWriter] Writing variable \'' + str(
                    VariableDefinition.LONG_NAME['longitude']) + '\' at time \'' + str(
                    time) + '\'')

            if self.points.is_raw_times:
                # Inutile de rechercher la date si l'axe du temps est le même que la donnée brute
                data[time_index] = self.points.read_variable_longitude_at_time(time_index)[self.index_x]
            else:
                data[time_index] = self.points.read_variable_longitude_at_time(time)[self.index_x]

            time_index += 1

        self.data['longitude'] = data

    def write_variable_latitude(self):
        logging.info('[DefaultTimePointWriter] Writing variable \''+str(VariableDefinition.LONG_NAME['latitude'])+'\'')
        data = np.zeros([self.points.get_t_size()])
        data[:] = np.nan
        time_index = 0
        for time in self.points.read_axis_t():

            logging.debug(
                '[DefaultTimePointWriter] Writing variable \'' + str(
                    VariableDefinition.LONG_NAME['latitude']) + '\' at time \'' + str(
                    time) + '\'')

            if self.points.is_raw_times:
                # Inutile de rechercher la date si l'axe du temps est le même que la donnée brute
                data[time_index] = self.points.read_variable_latitude_at_time(time_index)[self.index_x]
            else:
                data[time_index] = self.points.read_variable_latitude_at_time(time)[self.index_x]

            time_index += 1

        self.data['latitude'] = data

    def write_variable_sea_surface_wave_significant_height(self):
        logging.info('[DefaultTimePointWriter] Writing variable \''+str(VariableDefinition.LONG_NAME['sea_surface_wave_significant_height'])+'\'')
        data = np.zeros([self.points.get_t_size()])
        data[:] = np.nan
        time_index = 0
        for time in self.points.read_axis_t():

            logging.debug(
                '[DefaultTimePointWriter] Writing variable \'' + str(
                    VariableDefinition.LONG_NAME['sea_surface_wave_significant_height']) + '\' at time \'' + str(
                    time) + '\'')

            if self.points.is_raw_times:
                # Inutile de rechercher la date si l'axe du temps est le même que la donnée brute
                data[time_index] = self.points.read_variable_sea_surface_wave_significant_height_at_time(time_index)[
                    self.index_x]
            else:
                data[time_index] = self.points.read_variable_sea_surface_wave_significant_height_at_time(time)[self.index_x]

            time_index += 1

        self.data['sea_surface_wave_significant_height'] = data

    def write_variable_sea_surface_wave_mean_period(self):
        logging.info('[DefaultTimePointWriter] Writing variable \''+str(VariableDefinition.LONG_NAME['sea_surface_wave_mean_period'])+'\'')
        data = np.zeros([self.points.get_t_size()])
        data[:] = np.nan
        time_index = 0
        for time in self.points.read_axis_t():

            logging.debug(
                '[DefaultTimePointWriter] Writing variable \'' + str(
                    VariableDefinition.LONG_NAME['sea_surface_wave_mean_period']) + '\' at time \'' + str(
                    time) + '\'')

            if self.points.is_raw_times:
                # Inutile de rechercher la date si l'axe du temps est le même que la donnée brute
                data[time_index] = self.points.read_variable_sea_surface_wave_mean_period_at_time(time_index)[self.index_x]
            else:
                data[time_index] = self.points.read_variable_sea_surface_wave_mean_period_at_time(time)[self.index_x]

            time_index += 1

        self.data['sea_surface_wave_mean_periode'] = data

    def write_variable_sea_water_salinity_at_ground_level(self):
        logging.info('[DefaultTimePointWriter] Writing variable \''+str(VariableDefinition.LONG_NAME['sea_water_salinity_at_ground_level'])+'\'')
        data = np.zeros([self.points.get_t_size()])
        data[:] = np.nan
        time_index = 0
        for time in self.points.read_axis_t():

            logging.debug(
                '[DefaultTimePointWriter] Writing variable \'' + str(
                    VariableDefinition.LONG_NAME['sea_water_salinity_at_ground_level']) + '\' at time \'' + str(
                    time) + '\'')

            if self.points.is_raw_times:
                # Inutile de rechercher la date si l'axe du temps est le même que la donnée brute
                data[time_index] = self.points.read_variable_sea_water_salinity_at_ground_level_at_time(time_index)[
                    self.index_x]
            else:
                data[time_index] = self.points.read_variable_sea_water_salinity_at_ground_level_at_time(time)[self.index_x]

            time_index += 1

        self.data['sea_water_salinity_at_ground_level'] = data

    def write_variable_sea_surface_salinity(self):
        logging.info('[DefaultTimePointWriter] Writing variable \'' + str(
            VariableDefinition.LONG_NAME['sea_surface_salinity']) + '\'')
        data = np.zeros([self.points.get_t_size()])
        data[:] = np.nan
        time_index = 0
        for time in self.points.read_axis_t():

            logging.debug(
                '[DefaultTimePointWriter] Writing variable \'' + str(
                    VariableDefinition.LONG_NAME['sea_surface_salinity']) + '\' at time \'' + str(
                    time) + '\'')

            if self.points.is_raw_times:
                # Inutile de rechercher la date si l'axe du temps est le même que la donnée brute
                data[time_index] = self.points.read_variable_sea_surface_salinity_at_time(time_index)[self.index_x]
            else:
                data[time_index] = self.points.read_variable_sea_surface_salinity_at_time(time)[self.index_x]

            time_index += 1

        self.data['sea_surface_salinity'] = data

    def write_variable_sea_surface_density(self):
        logging.info('[DefaultTimePointWriter] Writing variable \''+str(VariableDefinition.LONG_NAME['sea_surface_density'])+'\'')
        data = np.zeros([self.points.get_t_size()])
        data[:] = np.nan
        time_index = 0
        for time in self.points.read_axis_t():

            logging.debug(
                '[DefaultTimePointWriter] Writing variable \'' + str(
                    VariableDefinition.LONG_NAME['sea_surface_density']) + '\' at time \'' + str(
                    time) + '\'')

            if self.points.is_raw_times:
                # Inutile de rechercher la date si l'axe du temps est le même que la donnée brute
                data[time_index] = self.points.read_variable_sea_surface_density_at_time(time_index)[self.index_x]
            else:
                data[time_index] = self.points.read_variable_sea_surface_density_at_time(time)[self.index_x]

            time_index += 1

        self.data['sea_surface_density'] = data

    def write_variable_sea_water_turbidity(self):
        logging.info('[DefaultTimePointWriter] Writing variable \''+str(VariableDefinition.LONG_NAME['sea_water_turbidity'])+'\'')
        data = np.zeros([self.points.get_t_size()])
        data[:] = np.nan
        time_index = 0
        for time in self.points.read_axis_t():

            logging.debug(
                '[DefaultTimePointWriter] Writing variable \'' + str(
                    VariableDefinition.LONG_NAME['sea_water_turbidity']) + '\' at time \'' + str(
                    time) + '\'')

            if self.points.is_raw_times:
                # Inutile de rechercher la date si l'axe du temps est le même que la donnée brute
                data[time_index] = self.points.read_variable_sea_water_turbidity_at_time(time_index)[self.index_x]
            else:
                data[time_index] = self.points.read_variable_sea_water_turbidity_at_time(time)[self.index_x]

            time_index += 1

        self.data['sea_water_turbidity'] = data

    def write_variable_bathymetry(self):
        logging.info('[DefaultTimePointWriter] Writing variable \''+str(VariableDefinition.LONG_NAME['bathymetry'])+'\'')
        data = np.zeros([self.points.get_t_size()])
        data[:] = np.nan
        time_index = 0
        for time in self.points.read_axis_t():

            logging.debug(
                '[DefaultTimePointWriter] Writing variable \'' + str(
                    VariableDefinition.LONG_NAME['bathymetry']) + '\' at time \'' + str(
                    time) + '\'')

            if self.points.is_raw_times:
                # Inutile de rechercher la date si l'axe du temps est le même que la donnée brute
                data[time_index] = self.points.read_variable_bathymetry_at_time(time_index)[self.index_x]
            else:
                data[time_index] = self.points.read_variable_bathymetry_at_time(time)[self.index_x]

            time_index += 1

        self.data['bathymetry'] = data

    def write_variable_wind_speed_10m(self):
        logging.info('[DefaultTimePointWriter] Writing variable \''+str(VariableDefinition.LONG_NAME['wind_speed_10m'])+'\'')
        data = np.zeros([self.points.get_t_size()])
        data[:] = np.nan
        time_index = 0
        for time in self.points.read_axis_t():

            logging.debug(
                '[DefaultTimePointWriter] Writing variable \'' + str(
                    VariableDefinition.LONG_NAME['wind_speed_10m']) + '\' at time \'' + str(
                    time) + '\'')

            if self.points.is_raw_times:
                # Inutile de rechercher la date si l'axe du temps est le même que la donnée brute
                data[time_index] = self.points.read_variable_wind_speed_10m_at_time(time_index)[self.index_x]
            else:
                data[time_index] = self.points.read_variable_wind_speed_10m_at_time(time)[self.index_x]

            time_index += 1

        self.data['wind_speed_10m'] = data

    def write_variable_wind_from_direction_10m(self):
        logging.info('[DefaultTimePointWriter] Writing variable \''+str(VariableDefinition.LONG_NAME['wind_from_direction_10m'])+'\'')
        data = np.zeros([self.points.get_t_size()])
        data[:] = np.nan
        time_index = 0
        for time in self.points.read_axis_t():

            logging.debug(
                '[DefaultTimePointWriter] Writing variable \'' + str(
                    VariableDefinition.LONG_NAME['wind_from_direction_10m']) + '\' at time \'' + str(
                    time) + '\'')

            if self.points.is_raw_times:
                # Inutile de rechercher la date si l'axe du temps est le même que la donnée brute
                data[time_index] = self.points.read_variable_wind_from_direction_10m_at_time(time_index)[self.index_x]
            else:
                data[time_index] = self.points.read_variable_wind_from_direction_10m_at_time(time)[self.index_x]

            time_index += 1

        self.data['wind_from_direction_10m'] = data

    def write_variable_wind_to_direction_10m(self):
        logging.info('[DefaultTimePointWriter] Writing variable \''+str(VariableDefinition.LONG_NAME['wind_to_direction_10m'])+'\'')
        data = np.zeros([self.points.get_t_size()])
        data[:] = np.nan
        time_index = 0
        for time in self.points.read_axis_t():

            logging.debug(
                '[DefaultTimePointWriter] Writing variable \'' + str(
                    VariableDefinition.LONG_NAME['wind_to_direction_10m']) + '\' at time \'' + str(
                    time) + '\'')

            if self.points.is_raw_times:
                # Inutile de rechercher la date si l'axe du temps est le même que la donnée brute
                data[time_index] = self.points.read_variable_wind_to_direction_10m_at_time(time_index)[self.index_x]
            else:
                data[time_index] = self.points.read_variable_wind_to_direction_10m_at_time(time)[self.index_x]

            time_index += 1

        self.data['wind_to_direction_10m'] = data

    def write_variable_surface_air_pressure(self):
        logging.info('[DefaultTimePointWriter] Writing variable \''+str(VariableDefinition.LONG_NAME['surface_air_pressure'])+'\'')
        data = np.zeros([self.points.get_t_size()])
        data[:] = np.nan
        time_index = 0
        for time in self.points.read_axis_t():

            logging.debug(
                '[DefaultTimePointWriter] Writing variable \'' + str(
                    VariableDefinition.LONG_NAME['surface_air_pressure']) + '\' at time \'' + str(
                    time) + '\'')

            if self.points.is_raw_times:
                # Inutile de rechercher la date si l'axe du temps est le même que la donnée brute
                data[time_index] = self.points.read_variable_surface_air_pressure_at_time(time_index)[self.index_x]
            else:
                data[time_index] = self.points.read_variable_surface_air_pressure_at_time(time)[self.index_x]

            time_index += 1

        self.data['surface_air_pressure'] = data

    def write_variable_rainfall_amount(self):
        logging.info('[DefaultTimePointWriter] Writing variable \''+str(VariableDefinition.LONG_NAME['rainfall_amount'])+'\'')
        data = np.zeros([self.points.get_t_size()])
        data[:] = np.nan
        time_index = 0
        for time in self.points.read_axis_t():

            logging.debug(
                '[DefaultTimePointWriter] Writing variable \'' + str(
                    VariableDefinition.LONG_NAME['rainfall_amount']) + '\' at time \'' + str(
                    time) + '\'')

            if self.points.is_raw_times:
                # Inutile de rechercher la date si l'axe du temps est le même que la donnée brute
                data[time_index] = self.points.read_variable_rainfall_amount_at_time(time_index)[self.index_x]
            else:
                data[time_index] = self.points.read_variable_rainfall_amount_at_time(time)[self.index_x]

            time_index += 1

        self.data['rainfall_amount'] = data

    def write_variable_sea_water_pressure_at_sea_water_surface(self):
        logging.info('[DefaultTimePointWriter] Writing variable \''+str(VariableDefinition.LONG_NAME['sea_water_pressure_at_sea_water_surface'])+'\'')
        data = np.zeros([self.points.get_t_size()])
        data[:] = np.nan
        time_index = 0
        for time in self.points.read_axis_t():

            logging.debug(
                '[DefaultTimePointWriter] Writing variable \'' + str(
                    VariableDefinition.LONG_NAME['sea_water_pressure_at_sea_water_surface']) + '\' at time \'' + str(
                    time) + '\'')

            if self.points.is_raw_times:
                # Inutile de rechercher la date si l'axe du temps est le même que la donnée brute
                data[time_index] = self.points.read_variable_sea_water_pressure_at_sea_water_surface_at_time(time_index)[
                    self.index_x]
            else:
                data[time_index] = self.points.read_variable_sea_water_pressure_at_sea_water_surface_at_time(time)[self.index_x]

            time_index += 1

        self.data['sea_water_pressure_at_sea_water_surface'] = data

    def write_variable_sea_water_electrical_conductivity(self):
        logging.info('[DefaultTimePointWriter] Writing variable \''+str(VariableDefinition.LONG_NAME['sea_water_electrical_conductivity'])+'\'')
        data = np.zeros([self.points.get_t_size()])
        data[:] = np.nan
        time_index = 0
        for time in self.points.read_axis_t():

            logging.debug(
                '[DefaultTimePointWriter] Writing variable \'' + str(
                    VariableDefinition.LONG_NAME['sea_water_electrical_conductivity']) + '\' at time \'' + str(
                    time) + '\'')

            if self.points.is_raw_times:
                # Inutile de rechercher la date si l'axe du temps est le même que la donnée brute
                data[time_index] = self.points.read_variable_sea_water_electrical_conductivity_at_time(time_index)[
                self.index_x]
            else:
                data[time_index] = self.points.read_variable_sea_water_electrical_conductivity_at_time(time)[
                self.index_x]

            time_index += 1

        self.data['sea_water_electrical_conductivity'] = data


    def write_variable_water_volume_transport_into_sea_water_from_rivers(self):
        logging.info('[DefaultTimePointWriter] Writing variable \''+str(VariableDefinition.LONG_NAME['water_volume_transport_into_sea_water_from_rivers'])+'\'')
        data = np.zeros([self.points.get_t_size()])
        data[:] = np.nan
        time_index = 0
        for time in self.points.read_axis_t():

            logging.debug(
                '[DefaultTimePointWriter] Writing variable \'' + str(
                    VariableDefinition.LONG_NAME['water_volume_transport_into_sea_water_from_rivers']) + '\' at time \'' + str(
                    time) + '\'')

            if self.points.is_raw_times:
                # Inutile de rechercher la date si l'axe du temps est le même que la donnée brute
                data[time_index] = self.points.read_variable_water_volume_transport_into_sea_water_from_rivers_at_time(time_index)[
                    self.index_x]
            else:
                data[time_index] = self.points.read_variable_water_volume_transport_into_sea_water_from_rivers_at_time(time)[
                self.index_x]

            time_index += 1

        self.data['water_volume_transport_into_sea_water_from_rivers'] = data
