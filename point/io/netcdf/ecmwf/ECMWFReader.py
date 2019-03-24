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
from point.io.MultiPointReader import MultiPointReader
from point.TimeMultiPoint import TimeMultiPoint
from coverage.TimeCoverage import TimeCoverage
from coverage.io.netcdf.ecmwf.ECMWFReader import ECMWFReader as CovReader
from netCDF4 import num2date
from scipy.io import loadmat
import numpy as np
from datetime import datetime,timedelta,timezone
import pytz
import logging


class ECMWFReader(MultiPointReader):

    def __init__(self,myFile,xy):
        MultiPointReader.__init__(self, myFile);
        self.coverage = TimeCoverage(CovReader(self.filename))
        self.nbPoints = np.shape(xy)[0]
        self.xy_coords = np.zeros([self.nbPoints,2],dtype=np.int32)
        self.xy_values = np.zeros([self.nbPoints, 2])
        self.meta_data = ""

        for i in range(0, self.nbPoints):
            nearestPoint = self.coverage.find_point_index(xy[i][0], xy[i][1])
            logging.info("Nearest point : " + str(nearestPoint[2]) + " / " + str(nearestPoint[3]) + " at " + str(round(nearestPoint[4], 4)) + " km")
            self.meta_data = self.meta_data + "\n# Point "+str(i)+" : nearest point in ECMWF file is "+ str(round(nearestPoint[4], 4)) + " km from the target point"
            logging.info("Nearest point (i,j) : " + str(nearestPoint[0]) + " / " + str(nearestPoint[1]))
            self.xy_coords[i] = [nearestPoint[0],nearestPoint[1]]
            self.xy_values[i] = [nearestPoint[2], nearestPoint[3]]

    def is_coverage_based(self):
            return True

    # Axis
    def read_axis_x(self):
        return self.xy_values[:,0]

    def read_axis_y(self):
        return self.xy_values[:,1]

    def read_axis_t(self,timestamp=0):
        return self.coverage.read_axis_t(timestamp)

    def read_metadata(self):
        m = {}
        m["data_source"] = "ECMWF file"
        m["meta_data"] = self.meta_data

        return m

    #Scalar
    def read_variable_wind_10m_at_time(self,date):

        data = np.zeros([self.nbPoints,2])
        data[:] = np.nan

        all_data = self.coverage.read_variable_wind_10m_at_time(date)

        for index_x in range(0, self.nbPoints):
            # comp U
            data[index_x,0] = all_data[0][self.xy_coords[index_x][1],self.xy_coords[index_x][0]]
            # comp V
            data[index_x,1] = all_data[1][self.xy_coords[index_x][1], self.xy_coords[index_x][0]]

        return data

    def read_variable_surface_air_pressure_at_time(self, date):

        data = np.zeros([self.nbPoints])
        data[:] = np.nan

        all_data = self.coverage.read_variable_surface_air_pressure_at_time(date)

        for index_x in range(0, self.nbPoints):
            data[index_x] = all_data[self.xy_coords[index_x][1], self.xy_coords[index_x][0]]

        return data

    def read_variable_rainfall_amount_at_time(self, date):

        data = np.zeros([self.nbPoints])
        data[:] = np.nan

        all_data = self.coverage.read_variable_rainfall_amount_at_time(date)

        for index_x in range(0, self.nbPoints):
            data[index_x] = all_data[self.xy_coords[index_x][1], self.xy_coords[index_x][0]]

        return data



