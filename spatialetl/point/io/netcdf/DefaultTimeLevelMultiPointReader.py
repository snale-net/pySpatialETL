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

import numpy as np
from netCDF4 import Dataset, num2date

from spatialetl.point.TimeMultiPoint import TimeMultiPoint
from spatialetl.point.io.MultiPointReader import MultiPointReader
from spatialetl.utils.VariableDefinition import VariableDefinition


class DefaultTimeLevelMultiPointReader(MultiPointReader):

    def __init__(self,myFile):
        MultiPointReader.__init__(self, myFile);
        self.ncfile = Dataset(self.filename, 'r')
        self.profilmax = np.shape(self.ncfile["point"])[0]
        self.zmax = np.shape(self.ncfile["depth"])[0]
        self.tmax = np.shape(self.ncfile["time"])[0]

    def get_z_size(self):
        return np.shape(self.ncfile["depth"])[0]

    # Axis
    def read_axis_x(self):
        return self.ncfile["longitude"]

    def read_axis_y(self):
        return self.ncfile["latitude"]

    def read_axis_z(self):
        return self.ncfile["depth"]

    def read_axis_t(self,timestamp=0):
        data = self.ncfile.variables['time'][:]
        result = num2date(data, units=self.ncfile.variables['time'].units, calendar="gregorian")

        if timestamp == 1:
            return [(t - TimeMultiPoint.TIME_DATUM).total_seconds() \
                    for t in result];
        else:
            return result;

    # Scalar
    def read_variable_sea_water_temperature_at_time_and_depth(self,index_t,index_z):

        return self.ncfile[VariableDefinition.VARIABLE_NAME['sea_water_temperature']][index_t,index_z]

    def read_variable_sea_water_salinity_at_time_and_depth(self,index_t,index_z):

        return self.ncfile[VariableDefinition.VARIABLE_NAME['sea_water_salinity']][index_t,index_z]

    def read_variable_sea_water_density_at_time_and_depth(self,index_t,index_z):

        return self.ncfile[VariableDefinition.VARIABLE_NAME['sea_water_density']][index_t,index_z]

    def read_variable_sea_water_turbidity_at_time_and_depth(self,index_t,index_z):

        return self.ncfile[VariableDefinition.VARIABLE_NAME['sea_water_turbidity']][index_t,index_z]

    def read_variable_sea_water_electrical_conductivity_at_time_and_depth(self,index_t,index_z):

        return self.ncfile[VariableDefinition.VARIABLE_NAME['sea_water_electrical_conductivity']][index_t,index_z]


