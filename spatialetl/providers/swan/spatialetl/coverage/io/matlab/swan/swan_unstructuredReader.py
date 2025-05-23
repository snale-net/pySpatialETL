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

import re

import cftime
import numpy as np
from scipy.io import loadmat

from spatialetl.coverage.io.CoverageReader import CoverageReader


class SWANUnstructuredReader(CoverageReader):

    def __init__(self,myFile):
        CoverageReader.__init__(self, myFile);
        self.mat = loadmat(self.filename)

        self.mat = loadmat(self.filename)

        temp = np.reshape(self.mat['Xp'][0], (-1, 1))
        self.x_size = np.shape(temp)[0]
        self.y_size = np.shape(temp)[1]

        self.times = []
        for var in self.mat:
            #print(var)
            if "Time" in var:
                t = re.search("Time_([0-9]{4})([0-9]{2})([0-9]{2})_([0-9]{2})([0-9]{2})([0-9]{2})", var)
                self.times.append(
                    cftime.datetime(int(t.group(1)), int(t.group(2)), int(t.group(3)), int(t.group(4)), int(t.group(5)),
                                    int(t.group(6))))

    def close(self):
        return

    def is_regular_grid(self):
        return False

    def get_x_size(self):
        return self.x_size

    def get_y_size(self):
        return self.y_size

    def get_t_size(self):
        return len(self.times)

    def read_axis_x(self,xmin,xmax,ymin,ymax):
        data = np.reshape(self.mat['Xp'][0], (self.y_size, self.x_size))[ymin:ymax,xmin:xmax]
        return data

    def read_axis_y(self, xmin, xmax, ymin, ymax):
        data = np.reshape(self.mat['Yp'][0], (self.y_size, self.x_size))[ymin:ymax, xmin:xmax]
        return data

    def read_variable_bathymetry(self, xmin, xmax, ymin, ymax):
        data = np.reshape(self.mat['Botlev'][0], (self.y_size, self.x_size))[ymin:ymax, xmin:xmax]
        return data

    def read_axis_t(self,tmin,tmax,timestamp):

        result = self.times[tmin:tmax]

        if timestamp == 1:
            return [(t - TimeCoverage.TIME_DATUM).total_seconds() \
                    for t in result];
        else:
            return result;

    def read_variable_sea_surface_wave_significant_height_at_time(self, t, xmin, xmax, ymin, ymax):
        time = self.times[t]
        key = 'Hsig_'+time.strftime("%Y%m%d_%H%M%S")

        data = np.reshape(self.mat[key][0], (self.y_size, self.x_size))[ymin:ymax, xmin:xmax]
        return data
        # try:
        #     if "hs" in self.ncfile.variables:
        #         return np.ma.filled(
        #             np.reshape(self.mat['Hsig'][0], (self.y_size, self.x_size))[ymin:ymax, xmin:xmax],
        #             fill_value=np.nan)
        # except Exception as ex:
        #     logging.debug("Error '" + str(ex) + "'")
        #     raise (VariableNameError("WW3UnstructuredReader", "An error occured : '" + str(ex) + "'", 1000))
        #
        # logging.debug("No variables found for '" + str(
        #     VariableDefinition.LONG_NAME['sea_surface_wave_significant_height']) + "'")
        # raise (VariableNameError("WW3UnstructuredReader",
        #                          "No variables found for '" + str(
        #                              VariableDefinition.LONG_NAME['sea_surface_wave_significant_height']) + "'",
        #                          1000))

    def read_variable_sea_surface_wave_mean_period_at_time(self, t, xmin, xmax, ymin, ymax):
        time = self.times[t]
        key = 'Tm01_'+time.strftime("%Y%m%d_%H%M%S")

        data = np.reshape(self.mat[key][0], (self.y_size, self.x_size))[ymin:ymax, xmin:xmax]
        return data
        # try:
        #     if "hs" in self.ncfile.variables:
        #         return np.ma.filled(
        #             np.reshape(self.mat['Hsig'][0], (self.y_size, self.x_size))[ymin:ymax, xmin:xmax],
        #             fill_value=np.nan)
        # except Exception as ex:
        #     logging.debug("Error '" + str(ex) + "'")
        #     raise (VariableNameError("WW3UnstructuredReader", "An error occured : '" + str(ex) + "'", 1000))
        #
        # logging.debug("No variables found for '" + str(
        #     VariableDefinition.LONG_NAME['sea_surface_wave_significant_height']) + "'")
        # raise (VariableNameError("WW3UnstructuredReader",
        #                          "No variables found for '" + str(
        #                              VariableDefinition.LONG_NAME['sea_surface_wave_significant_height']) + "'",
        #                          1000))

    def read_variable_sea_surface_wave_from_direction_at_time(self, t, xmin, xmax, ymin, ymax):
        time = self.times[t]
        key = 'Dir_' + time.strftime("%Y%m%d_%H%M%S")

        data = 270.-np.reshape(self.mat[key][0], (self.y_size, self.x_size))[ymin:ymax, xmin:xmax]
        return data
        # try:
        #     if "hs" in self.ncfile.variables:
        #         return np.ma.filled(
        #             np.reshape(self.mat['Hsig'][0], (self.y_size, self.x_size))[ymin:ymax, xmin:xmax],
        #             fill_value=np.nan)
        # except Exception as ex:
        #     logging.debug("Error '" + str(ex) + "'")
        #     raise (VariableNameError("WW3UnstructuredReader", "An error occured : '" + str(ex) + "'", 1000))
        #
        # logging.debug("No variables found for '" + str(
        #     VariableDefinition.LONG_NAME['sea_surface_wave_significant_height']) + "'")
        # raise (VariableNameError("WW3UnstructuredReader",
        #                          "No variables found for '" + str(
        #                              VariableDefinition.LONG_NAME['sea_surface_wave_significant_height']) + "'",
        #                          1000))

    def read_variable_sea_surface_wave_to_direction_at_time(self, t, xmin, xmax, ymin, ymax):
        time = self.times[t]
        key = 'Dir_' + time.strftime("%Y%m%d_%H%M%S")

        data = np.reshape(self.mat[key][0], (self.y_size, self.x_size))[ymin:ymax, xmin:xmax]
        return data
        # try:
        #     if "hs" in self.ncfile.variables:
        #         return np.ma.filled(
        #             np.reshape(self.mat['Hsig'][0], (self.y_size, self.x_size))[ymin:ymax, xmin:xmax],
        #             fill_value=np.nan)
        # except Exception as ex:
        #     logging.debug("Error '" + str(ex) + "'")
        #     raise (VariableNameError("WW3UnstructuredReader", "An error occured : '" + str(ex) + "'", 1000))
        #
        # logging.debug("No variables found for '" + str(
        #     VariableDefinition.LONG_NAME['sea_surface_wave_significant_height']) + "'")
        # raise (VariableNameError("WW3UnstructuredReader",
        #                          "No variables found for '" + str(
        #                              VariableDefinition.LONG_NAME['sea_surface_wave_significant_height']) + "'",
        #                          1000))

    def read_variable_wind_10m_at_time(self,t,xmin,xmax,ymin,ymax):
        time = self.times[t]
        u_key = 'Windv_x_' + time.strftime("%Y%m%d_%H%M%S")
        v_key = 'Windv_y_' + time.strftime("%Y%m%d_%H%M%S")
        return [np.reshape(self.mat[u_key][0], (self.y_size, self.x_size))[ymin:ymax,xmin:xmax],
                np.reshape(self.mat[v_key][0], (self.y_size, self.x_size))[ymin:ymax,xmin:xmax] ]

