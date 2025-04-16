#! /usr/bin/env python2.7
# -*- coding: utf-8 -*-
# MIT License
# Copyright (c) 2024 [SNALE - French SAS Company - RCS 951 724 616]
#
# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:
#
# The above copyright notice and this permission notice shall be included in all
# copies or substantial portions of the Software.
#
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
# SOFTWARE.
from __future__ import division, print_function, absolute_import

import math
import re
import numpy as np

import pandas

from spatialetl.point.TimeMultiPoint import TimeMultiPoint
from spatialetl.point.io.MultiPointReader import MultiPointReader
from spatialetl.utils.logger import logging


class DefaultTimePointReader(MultiPointReader):
    def __init__(self, myFilename,colsNumber,varNames,names=[],checkOverlapping=False):
        """
        Initialize a reader
        :param myFilename: Path of the file
        :param colsNumber: List of index cols to read
        :param varNames: List of parameter names in the same order than colsNumber
        :param names: (Optional) List of point names
        :param checkOverlapping: (Optional) Check if overlapping
        """

        if isinstance(myFilename,str):
            MultiPointReader.__init__(self,myFilename)
        elif isinstance(myFilename,list):
            MultiPointReader.__init__(self, myFilename[0])

        self.header = 0
        self.x = ["Undefinied"]
        self.y = ["Undefinied"]
        self.names = ["Undefinied"]
        self.read_metadata()

        if len(names) != 0 :
            self.names = names

        if isinstance(myFilename, str):

            self.data = pandas.read_csv(self.filename, usecols=colsNumber,
                                   names=varNames, sep='\t', na_values={"NaN"},
                                   keep_default_na=False, parse_dates={
                    'time': ['date']},header=self.header)

            self.data = self.data.set_index(pandas.DatetimeIndex(self.data['time']))
            self.data = self.data.drop('time', 1)

        else:
            li = []
            for candidateFile in myFilename:

                df = pandas.read_csv(candidateFile, usecols=colsNumber,
                                            names=varNames,
                                            sep='\t', na_values={"NaN"},
                                            keep_default_na=False, parse_dates={
                        'time': ['date']}, header=self.header)
                df = df.set_index(pandas.DatetimeIndex(df['time']))
                df = df.drop('time', 1)

                li.append(df)

            try:
                self.data = pandas.concat(li, axis=0,verify_integrity=checkOverlapping)
            except ValueError as ex:

                if "overlapping" in str(ex.args):
                    logging.info("[DefaultTimePointReader] Dates from multiple files "+str(myFilename)+" are overlapping. Abort.")
                    raise ValueError("[DefaultTimePointReader] Dates from multiple files are overlapping. Abort.")
                else:
                    raise ValueError(ex)

    # Axis
    def get_x_size(self):
        return len(self.x);

    def get_y_size(self):
        return len(self.y);

    def get_t_size(self):
        return len(self.data.index.to_pydatetime());

    def read_axis_x(self):
        return self.x

    def read_axis_y(self):
        return self.y

    def read_axis_t(self,tmin,tmax,timestamp=0):
        """
        Returns the time axis
        :param tmin: Min slice index
        :param tmax: Max slice index
        :param timestamp: True if time in timestamp, else datime
        :return: Array of time
        """

        result = self.data.index.to_pydatetime()

        if timestamp == 1:
            return [(t - TimeMultiPoint.TIME_DATUM).total_seconds() \
                    for t in result][tmin:tmax];
        else:
            return result[tmin:tmax]

    def read_metadata(self):
        metadata = {};
        self.header = -1
        with open(self.filename) as f:
            line = f.readline()
            while line and line.startswith("#"):

                self.header = self.header +1

                # TODO : gérer les espaces blancs

                if "Station" in line:
                    metadata['name_station'] = re.sub('[^a-zA-Z0-9-_*.]', '', line.rsplit(':', 1)[1])
                    self.names = [metadata['name_station']]

                if "Longitude" in line:
                    metadata['x_coord'] = re.sub('[^a-zA-Z0-9-_*.]', '', line.rsplit(':', 1)[1])
                    self.x = [float(metadata['x_coord'])]

                if "Latitude" in line:
                    metadata['y_coord'] = re.sub('[^a-zA-Z0-9-_*.]', '', line.rsplit(':', 1)[1])
                    self.y = [float(metadata['y_coord'])]

                if "Vertical datum" in line:
                    metadata['vertical_datum'] = re.sub('[^a-zA-Z0-9-_*.]', '', line.rsplit(':', 1)[1])

                if "Data source" in line:
                    metadata['data_source'] = re.sub('[^a-zA-Z0-9-_*.]', '', line.rsplit(':', 1)[1])

                line = f.readline()

        return metadata

    def read_variable_point_names(self):
        return self.names

    #################
    # HYDRO
    # Sea Surface
    #################
    def read_variable_sea_surface_height_above_mean_sea_level_at_time(self, index_t):

        result = self.data.iloc[index_t].sea_surface_height_above_mean_sea_level

        return [result]

    def read_variable_sea_water_column_thickness_at_time(self, index_t):

        result = self.data.iloc[index_t].sea_water_column_thickness

        return [result]

    def read_variable_sea_water_pressure_at_sea_water_surface_at_time(self, index_t):

        result = self.data.iloc[index_t].sea_water_pressure_at_sea_water_surface

        return [result]

    #################
    # HYDRO
    # Ground level
    #################

    #################
    # HYDRO
    # 2D
    #################
    def read_variable_water_volume_transport_into_sea_water_from_rivers_at_time(self, index_t):

        result = self.data.iloc[index_t].water_volume_transport_into_sea_water_from_rivers

        return [result]

    #################
    # WAVES
    # Sea Surface
    #################
    def read_variable_sea_surface_wave_significant_height_at_time(self, index_t):

        result = self.data.iloc[index_t].sea_surface_wave_significant_height

        return [result]

    #################
    # WAVES
    # Momentum flux
    #################

    #################
    # METEO
    # 2D
    #################

    #################
    # METEO
    # Sea surface
    #################
    def read_variable_surface_air_pressure_at_time(self, index_t):

        result = self.data.iloc[index_t].surface_air_pressure

        return [result]

    def read_variable_sea_surface_air_pressure_at_time(self, index_t):

        result = self.data.iloc[index_t].sea_surface_air_pressure

        return [result]

    def read_variable_rainfall_amount_at_time(self, index_t):

        result = self.data.iloc[index_t].rainfall_amount
        # to kg/m²
        # result = self.data.iloc[index_t].rainfall_amount

        return [result]

    #################
    # METEO
    # At 10 m
    #################
    def read_variable_wind_10m_at_time(self,index_t):

        u = self.data.iloc[index_t].wind_speed_10m * math.cos(math.radians(self.data.iloc[index_t].wind_from_direction_10m))
        v = self.data.iloc[index_t].wind_speed_10m * math.sin(math.radians(self.data.iloc[index_t].wind_from_direction_10m))

        return [[u],[v]]
