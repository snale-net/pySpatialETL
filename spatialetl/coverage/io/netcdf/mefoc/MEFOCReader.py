#! /usr/bin/env python3.6
# -*- coding: utf-8 -*-
# Author : Fabien Rétif - fabien.retif@zoho.com
#
from __future__ import division, print_function, absolute_import

import os
from datetime import datetime

import numpy as np
from netCDF4 import Dataset, MFDataset, num2date

from spatialetl.coverage.TimeCoverage import TimeCoverage
from spatialetl.coverage.io.CoverageReader import CoverageReader


class MEFOCReader(CoverageReader):

    def __init__(self, myFile):
        CoverageReader.__init__(self, myFile);

        if os.path.isfile(self.filename):
            self.ncfile = Dataset(self.filename, 'r')
        elif os.path.isdir(self.filename):
            self.ncfile = MFDataset(os.path.join(self.filename, "*.nc"), 'r')
        elif self.filename.endswith("*"):
            self.ncfile = MFDataset(self.filename + ".nc", 'r')
        else:
            raise ValueError("Unable to decode file " + str(self.filename))

    def close(self):
        self.ncfile.close()

    def is_regular_grid(self):
        return True

    def get_x_size(self):
        return np.shape(self.ncfile.variables['longitude'][:])[0];

    def get_y_size(self):
        return np.shape(self.ncfile.variables['latitude'][:])[0];

    def get_z_size(self):
        return np.shape(self.ncfile.variables['depth'][:])[0];

    def get_t_size(self):
        return np.shape(self.ncfile.variables['time'])[0];

    # Axis
    def read_axis_t(self,tmin,tmax, timestamp=0):
        data = self.ncfile.variables['time'][tmin:tmax]
        temp = num2date(data, units=self.ncfile.variables['time'].units, calendar=self.ncfile.variables['time'].calendar)

        result = [datetime.strptime(str(t), '%Y-%m-%d %H:%M:%S') \
                  for t in temp];

        if timestamp == 1:
            return [(t - TimeCoverage.TIME_DATUM).total_seconds() \
                    for t in result];
        else:
            return result;

    def read_axis_x(self,xmin,xmax,ymin,ymax):
        return self.ncfile.variables['longitude'][xmin:xmax]

    def read_axis_y(self,xmin,xmax,ymin,ymax):
        return self.ncfile.variables['latitude'][ymin:ymax]

    # Scalar
    def read_variable_bathymetry(self,xmin,xmax,ymin,ymax):
        return np.ma.filled(self.ncfile.variables['bathymetry'][ymin:ymax,xmin:xmax], fill_value=np.nan)

    def read_variable_eastward_wind_10m(self,xmin,xmax,ymin,ymax):
        return np.ma.filled(self.ncfile.variables["u_wind_10m"][:,ymin:ymax,xmin:xmax], fill_value=np.nan)

    def read_variable_northward_wind_10m(self,xmin,xmax,ymin,ymax):
        return np.ma.filled(self.ncfile.variables["v_wind_10m"][:,ymin:ymax,xmin:xmax], fill_value=np.nan)

    def read_variable_wind_10m_at_time(self, index_t,xmin,xmax,ymin,ymax):
        return [np.ma.filled(self.ncfile.variables["u_wind_10m"][index_t][ymin:ymax,xmin:xmax], fill_value=np.nan),
                np.ma.filled(self.ncfile.variables["v_wind_10m"][index_t][ymin:ymax,xmin:xmax], fill_value=np.nan)]

    def read_variable_sea_surface_wave_significant_height(self,xmin,xmax,ymin,ymax):
            return self.ncfile.variables["hs"][::]

    def read_variable_sea_surface_wave_mean_period(self,xmin,xmax,ymin,ymax):
            return self.ncfile.variables["tm0m1"][::]

    def read_variable_sea_surface_wave_to_direction(self):
        return self.ncfile.variables["dir"][::]

