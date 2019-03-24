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
from coverage.TimeLevelCoverage import TimeLevelCoverage
from coverage.io.netcdf.symphonie.SymphonieReader import SymphonieReader as CovReader
from netCDF4 import num2date
from scipy.io import loadmat
import numpy as np
from datetime import datetime,timedelta,timezone
import pytz
import logging


class SymphonieReader(MultiPointReader):

    def __init__(self,myGrid,myFile,xy,depth):
        MultiPointReader.__init__(self, myFile);
        self.coverage = TimeLevelCoverage(CovReader(myGrid,self.filename))
        self.nbPoints = np.shape(xy)[0]
        self.xy_coords = np.zeros([self.nbPoints,2],dtype=np.int32)
        self.xy_values = np.zeros([self.nbPoints, 2])
        self.meta_data = ""

        for i in range(0, np.shape(self.xy_coords)[0]):
            nearestPoint = self.coverage.find_point_index(xy[i][0], xy[i][1])
            logging.info("Nearest point : " + str(nearestPoint[2]) + " / " + str(nearestPoint[3]) + " at " + str(round(nearestPoint[4], 4)) + " km")
            self.meta_data = self.meta_data + "\n# Point " + str(
                i) + " : nearest point in SYMPHONIE file is " + str(round(nearestPoint[4],4)) + " km from the target point"
            logging.info("Nearest point (i,j) : " + str(nearestPoint[0]) + " / " + str(nearestPoint[1]))
            self.xy_coords[i] = [nearestPoint[0],nearestPoint[1]]
            self.xy_values[i] = [nearestPoint[2], nearestPoint[3]]

        self.depths = depth

    def is_coverage_based(self):
        return True

    # Axis
    def read_axis_x(self):
        return self.xy_values[:,0]

    def read_axis_y(self):
        return self.xy_values[:,1]

    def read_axis_z(self):
        return self.depths

    def read_axis_t(self,timestamp=0):
        return self.coverage.read_axis_t(timestamp)

    def read_metadata(self):
        m = {}
        m["data_source"] = "SYMPHONIE file"
        m["meta_data"] = self.meta_data

        return m

    #Scalar
    def read_variable_sea_water_temperature_at_time_and_depth(self,index_t,depth):

        data = np.zeros([self.nbPoints])
        data[:] = np.nan

        all_data = self.coverage.read_variable_sea_water_temperature_at_time_and_depth(index_t,depth)

        for index_x in range(0, self.nbPoints):
            data[index_x] = all_data[self.xy_coords[index_x][1],self.xy_coords[index_x][0]]

        return data

    def read_variable_sea_water_salinity_at_time_and_depth(self, index_t, depth):

        data = np.zeros([self.nbPoints])
        data[:] = np.nan

        all_data = self.coverage.read_variable_sea_water_salinity_at_time_and_depth(index_t, depth)

        for index_x in range(0, self.nbPoints):
            data[index_x] = all_data[self.xy_coords[index_x][1], self.xy_coords[index_x][0]]

        return data

    def read_variable_baroclinic_sea_water_velocity_at_time_and_depth(self, index_t, depth):

        data = np.zeros([2,self.nbPoints])
        data[:] = np.nan

        all_data = self.coverage.read_variable_baroclinic_sea_water_velocity_at_time_and_depth(index_t, depth)

        for index_x in range(0, self.nbPoints):
            data[0][index_x] = all_data[0][self.xy_coords[index_x][1], self.xy_coords[index_x][0]]
            data[1][index_x] = all_data[1][self.xy_coords[index_x][1], self.xy_coords[index_x][0]]

        return data

    def read_variable_barotropic_sea_water_velocity_at_time(self, index_t):

        data = np.zeros([2,self.nbPoints])
        data[:] = np.nan

        all_data = self.coverage.read_variable_barotropic_sea_water_velocity_at_time(index_t)

        for index_x in range(0, self.nbPoints):
            data[0][index_x] = all_data[0][self.xy_coords[index_x][1], self.xy_coords[index_x][0]]
            data[1][index_x] = all_data[1][self.xy_coords[index_x][1], self.xy_coords[index_x][0]]

        return data

    def read_variable_sea_surface_temperature_at_time(self, index_t):

        depth = 0.0
        data = np.zeros([self.nbPoints])
        data[:] = np.nan

        all_data = self.coverage.read_variable_sea_water_temperature_at_time_and_depth(index_t, depth)

        for index_x in range(0, self.nbPoints):
            data[index_x] = all_data[self.xy_coords[index_x][1], self.xy_coords[index_x][0]]

        return data

    def read_variable_sea_surface_salinity_at_time(self, index_t):

        depth = 0.0
        data = np.zeros([self.nbPoints])
        data[:] = np.nan

        all_data = self.coverage.read_variable_sea_water_salinity_at_time_and_depth(index_t, depth)

        for index_x in range(0, self.nbPoints):
            data[index_x] = all_data[self.xy_coords[index_x][1], self.xy_coords[index_x][0]]

        return data




