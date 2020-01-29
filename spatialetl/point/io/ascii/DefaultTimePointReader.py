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
from spatialetl.point.io.MultiPointReader import MultiPointReader
from spatialetl.point.TimeMultiPoint import TimeMultiPoint
import os.path
import pandas
import logging
import math
import numpy as np
from datetime import datetime
import pytz
import re

class DefaultTimePointReader(MultiPointReader):
    def __init__(self, myFilename,names,colsNumber,varNames,checkOverlapping=False):

        if isinstance(myFilename,str):
            MultiPointReader.__init__(self,myFilename)
        elif isinstance(myFilename,list):
            MultiPointReader.__init__(self, myFilename[0])

        self.names = names
        self.header = 0
        self.x = ["Undefinied"]
        self.y = ["Undefinied"]
        self.read_metadata()

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
    def read_axis_x(self):
        return self.x

    def read_axis_y(self):
        return self.y

    def read_axis_t(self,timestamp=0):

        result = self.data.index.to_pydatetime()

        if timestamp == 1:
            return [(t - TimeMultiPoint.TIME_DATUM).total_seconds() \
                    for t in result];
        else:
            return result

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

                if "Longitude" in line:
                    metadata['x_coord'] = re.sub('[^a-zA-Z0-9-_*.]', '', line.rsplit(':', 1)[1])
                    self.x = [metadata['x_coord']]

                if "Latitude" in line:
                    metadata['y_coord'] = re.sub('[^a-zA-Z0-9-_*.]', '', line.rsplit(':', 1)[1])
                    self.y = [metadata['y_coord']]

                if "Vertical datum" in line:
                    metadata['vertical_datum'] = re.sub('[^a-zA-Z0-9-_*.]', '', line.rsplit(':', 1)[1])

                if "Data source" in line:
                    metadata['data_source'] = re.sub('[^a-zA-Z0-9-_*.]', '', line.rsplit(':', 1)[1])

                line = f.readline()

        return metadata

    def read_variable_wind_10m_at_time(self,index_t):

        u = self.data.iloc[index_t].wind_speed_10m * math.cos(math.radians(self.dir_data.iloc[index_t].wind_from_direction_10m))
        v = self.data.iloc[index_t].wind_speed_10m * math.sin(math.radians(self.dir_data.iloc[index_t].wind_from_direction_10m))

        return [[u],[v]]

    def read_variable_surface_air_pressure_at_time(self,index_t):

        # bar to Pa
        result = self.data.iloc[index_t].surface_air_pressure * 100000

        return [result]

    def read_variable_rainfall_amount_at_time(self,index_t):

        result = self.data.iloc[index_t].rainfall_amount
        # to kg/m²
        #result = self.data.iloc[index_t].rainfall_amount

        return [result]

    def read_variable_point_names(self):
        return self.names

    def read_variable_water_volume_transport_into_sea_water_from_rivers_at_time(self,index_t):

        result = self.data.iloc[index_t].water_volume_transport_into_sea_water_from_rivers

        return [result]

    def read_variable_sea_water_pressure_at_sea_water_surface_at_time(self, index_t):

        result = self.data.iloc[index_t].sea_water_pressure_at_sea_water_surface

        return [result]

    def read_variable_sea_surface_height_above_mean_sea_level_at_time(self, index_t):

        result = self.data.iloc[index_t].sea_surface_height_above_mean_sea_level

        return [result]
