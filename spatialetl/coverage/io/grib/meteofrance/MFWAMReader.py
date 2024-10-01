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

from datetime import datetime
from datetime import timedelta

import cfgrib
import numpy as np

from spatialetl.coverage.TimeCoverage import TimeCoverage
from spatialetl.coverage.io.CoverageReader import CoverageReader
from spatialetl.exception.VariableNameError import VariableNameError
from spatialetl.utils.VariableDefinition import VariableDefinition
from spatialetl.utils.logger import logging


class MFWAMReader (CoverageReader):

    def __init__(self, myFile):
        CoverageReader.__init__(self,myFile);

        try:
            self.file = cfgrib.open_file(self.filename)
        except Exception as ex:
            logging.debug("Error '" + str(ex) + "'")
            raise VariableNameError("MFWAMReader", "An error occured : '" + str(ex) + "'", 1000)

        logging.debug(sorted(self.file.variables))

        self.t_size = 1
        self.times = [self.file.variables['time'].data]

        lon = self.file.variables['longitude'].data
        self.new_lon = np.mod(lon + 180.0, 360.0) - 180.0
        lat = self.file.variables['latitude'].data
        xx, yy = np.meshgrid(self.new_lon, lat)
        self.new_indexes = np.argsort(xx)

    def find_time_and_step(self,index_t):
        try:
            if "time" in self.file.variables and "step" in self.file.variables:
                return 0, 0
                array = np.asarray(self.file.variables["time"].data)
                idx = np.where(array <= self.file.variables["time"].data[int(index_t)])

                if len(idx[0]) == 0:
                    raise ValueError("Time "+self.file.variables["time"].data[int(index_t)]+" was not found")

                nearest_t_index = (np.abs(array[idx] - self.file.variables["time"].data[int(index_t)])).argmin()

                td = datetime.utcfromtimestamp(
                    self.file.variables["time"].data[int(index_t)]) - datetime.utcfromtimestamp(
                    self.file.variables["time"].data[int(nearest_t_index)])

                array = np.asarray(self.ds_accum.variables["step"].data)
                nearest_step_index = (np.abs(array - divmod(td.seconds, 3600)[0])).argmin()

                #logging.debug("On cherche :",datetime.utcfromtimestamp(self.ds_instep.variables["time"].data[int(index_t)]) )
                #logging.debug("ref Time",datetime.utcfromtimestamp(self.ds_step.variables["time"].data[int(nearest_t_index)]))
                step = timedelta(hours=self.file.variables["step"].data[int(nearest_step_index)])
                #logging.debug("On trouve",datetime.utcfromtimestamp(self.ds_step.variables["time"].data[int(nearest_t_index)])+step)

                return nearest_t_index, nearest_step_index
            else:
                logging.debug("No variables found for '" + str(VariableDefinition.LONG_NAME['time']) + "'")
                raise (VariableNameError("ECMWFReader",
                                         "No variables found for '" + str(
                                             VariableDefinition.LONG_NAME['time']) + "'",
                                         1000))
        except Exception as ex:
            logging.debug("Error '" + str(ex) + "'")
            raise (VariableNameError("ECMWFReader", "An error occured : '" + str(ex) + "'", 1000))

    def close(self):
        self.file.close()

    def is_regular_grid(self):
        return True

    def get_x_size(self):
        try:
            if "longitude" in self.file.variables:
                return self.file.dimensions["longitude"]
            else:
                logging.debug("No variables found for '" + str(VariableDefinition.LONG_NAME['longitude']) + "'")
                raise (VariableNameError("ECMWFReader",
                                         "No variables found for '" + str(
                                             VariableDefinition.LONG_NAME['longitude']) + "'",
                                         1000))
        except Exception as ex:
            logging.debug("Error '" + str(ex) + "'")
            raise (VariableNameError("ECMWFReader", "An error occured : '" + str(ex) + "'", 1000))

    def get_y_size(self):
        try:
            if "latitude" in self.file.variables:
                return self.file.dimensions["latitude"]
            else:
                logging.debug("No variables found for '" + str(VariableDefinition.LONG_NAME['latitude']) + "'")
                raise (VariableNameError("ECMWFReader",
                                         "No variables found for '" + str(
                                             VariableDefinition.LONG_NAME['latitude']) + "'",
                                         1000))
        except Exception as ex:
            logging.debug("Error '" + str(ex) + "'")
            raise (VariableNameError("ECMWFReader", "An error occured : '" + str(ex) + "'", 1000))

    def get_t_size(self):
        return self.t_size

    def read_axis_x(self,xmin,xmax,ymin,ymax):
        return np.array(sorted(self.new_lon)[xmin:xmax])

    def read_axis_y(self,xmin,xmax,ymin,ymax):
        try:
            if "latitude" in self.file.variables:
                return self.file.variables['latitude'].data[ymin:ymax]
            else:
                logging.debug("No variables found for '" + str(VariableDefinition.LONG_NAME['latitude']) + "'")
                raise (VariableNameError("ECMWFReader",
                                         "No variables found for '" + str(
                                             VariableDefinition.LONG_NAME['latitude']) + "'",
                                         1000))
        except Exception as ex:
            logging.debug("Error '" + str(ex) + "'")
            raise (VariableNameError("ECMWFReader", "An error occured : '" + str(ex) + "'", 1000))

    def read_axis_t(self, tmin, tmax, timestamp):

        data = self.times[tmin:tmax]
        result = [datetime.utcfromtimestamp(t) \
                  for t in data]

        if timestamp == 1:
            return [(t - TimeCoverage.TIME_DATUM).total_seconds() \
                    for t in result];
        else:
            return result;

    # Variables
    def read_variable_longitude(self,xmin,xmax,ymin,ymax):
        return self.read_axis_x()

    def read_variable_latitude(self,xmin,xmax,ymin,ymax):
        return self.read_axis_y()

    def read_variable_time(self,tmin,tmax):
        return self.read_axis_t(tmin,tmax,timestamp=0)

    #################
    # WAVES
    # Sea Surface
    #################

    def read_variable_sea_surface_wave_significant_height_at_time(self, index_t,xmin,xmax,ymin,ymax):
        try:
            if "p3100" in self.file.variables:
                nearest_t_index, nearest_step_index = self.find_time_and_step(index_t)
                data = np.ma.filled(self.file.variables["p3100"].data[int(nearest_t_index), int(nearest_step_index),
                                    int(ymin):int(ymax), :], fill_value=np.nan)
                return np.take_along_axis(data, self.new_indexes[ymin:ymax], axis=1)[:, xmin:xmax]
            else:
                logging.debug(
                    "No variables found for '" + str(VariableDefinition.LONG_NAME['rainfall_amount']) + "'")
                raise (VariableNameError("MFWAMReader",
                                         "No variables found for '" + str(
                                             VariableDefinition.LONG_NAME['rainfall_amount']) + "'",
                                         1000))
        except Exception as ex:
            logging.debug("Error '" + str(ex) + "'")
            raise (VariableNameError("MFWAMReader", "An error occured : '" + str(ex) + "'", 1000))




