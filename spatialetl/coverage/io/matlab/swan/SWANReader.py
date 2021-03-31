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

import re

import cftime
import numpy as np
from scipy.io import loadmat

from spatialetl.coverage.TimeCoverage import TimeCoverage
from spatialetl.coverage.io.CoverageReader import CoverageReader


class SWANReader(CoverageReader):

    def __init__(self,myFile):
        CoverageReader.__init__(self, myFile);
        self.mat = loadmat(self.filename)

        self.mat = loadmat(self.filename)

        self.x_size = np.shape(self.mat['Xp'])[1]
        self.y_size = np.shape(self.mat['Yp'])[0]
        self.times = []
        for var in self.mat:
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
        data = self.mat['Xp'][ymin:ymax,xmin:xmax]
        return data

    def read_axis_y(self, xmin, xmax, ymin, ymax):
        data = self.mat['Yp'][ymin:ymax, xmin:xmax]
        return data

    def read_variable_bathymetry(self, xmin, xmax, ymin, ymax):
        data = self.mat['Botlev'][ymin:ymax, xmin:xmax]
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

        data = self.mat[key][ymin:ymax, xmin:xmax]
        return data

    def read_variable_sea_surface_wave_mean_period_at_time(self, t, xmin, xmax, ymin, ymax):
        time = self.times[t]
        key = 'Tm01_'+time.strftime("%Y%m%d_%H%M%S")

        data = self.mat[key][ymin:ymax, xmin:xmax]
        return data

    def read_variable_sea_surface_wave_from_direction_at_time(self, t, xmin, xmax, ymin, ymax):
        time = self.times[t]
        key = 'Dir_' + time.strftime("%Y%m%d_%H%M%S")

        data = 270.-self.mat[key][ymin:ymax, xmin:xmax]
        return data

    def read_variable_sea_surface_wave_to_direction_at_time(self, t, xmin, xmax, ymin, ymax):
        time = self.times[t]
        key = 'Dir_' + time.strftime("%Y%m%d_%H%M%S")

        data = self.mat[key][ymin:ymax, xmin:xmax]
        return data

    def read_variable_wind_10m_at_time(self,t,xmin,xmax,ymin,ymax):
        time = self.times[t]
        u_key = 'Windv_x_' + time.strftime("%Y%m%d_%H%M%S")
        v_key = 'Windv_y_' + time.strftime("%Y%m%d_%H%M%S")
        return [self.mat[u_key][ymin:ymax,xmin:xmax],
                self.mat[v_key][ymin:ymax,xmin:xmax]]

