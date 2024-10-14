#! /usr/bin/env python3.6
# -*- coding: utf-8 -*-
# Author : Fabien Rétif - fabien.retif@snale.net
#
from __future__ import division, print_function, absolute_import

import os
from datetime import datetime

import numpy as np
from netCDF4 import Dataset, MFDataset, num2date

from spatialetl.coverage.TimeCoverage import TimeCoverage
from spatialetl.coverage.io.CoverageReader import CoverageReader


class HYCOMReader(CoverageReader):

    def __init__(self,myGrid, myFile=None):
        CoverageReader.__init__(self, myGrid);

        self.grid = Dataset(self.filename, 'r')

        if os.path.isfile(myFile):
            self.ncfile = Dataset(myFile, 'r')
        elif os.path.isdir(myFile):
            self.ncfile = MFDataset(os.path.join(myFile, "*.nc"), 'r')
        elif myFile.endswith("*"):
            self.ncfile = MFDataset(myFile + ".nc", 'r')
        else:
            raise ValueError("Unable to decode file " + str(myFile))

    def close(self):
        self.ncfile.close()

    def is_regular_grid(self):
        return True

    def get_x_size(self):
        return np.shape(self.grid.variables['lon'][:])[0];

    def get_y_size(self):
        return np.shape(self.grid.variables['lat'][:])[0];

    def get_t_size(self):
        return np.shape(self.ncfile.variables['time'])[0];

    # Axis
    def read_axis_t(self,tmin,tmax, timestamp=0):
        data = self.ncfile.variables['time'][tmin:tmax]
        temp = num2date(data, units=self.ncfile.variables['time'].units, calendar='gregorian')

        result = [datetime.strptime(str(t), '%Y-%m-%d %H:%M:%S') \
                  for t in temp];

        if timestamp == 1:
            return [(t - TimeCoverage.TIME_DATUM).total_seconds() \
                    for t in result];
        else:
            return result;

    def read_axis_x(self,xmin,xmax,ymin,ymax):
        return self.grid.variables['lon'][ymin:ymax, xmin:xmax]

    def read_axis_y(self,xmin,xmax,ymin,ymax):
        return self.grid.variables['lat'][ymin:ymax, xmin:xmax]

    #################
    # HYDRO
    # 2D
    #################
    def read_variable_bathymetry(self,xmin,xmax,ymin,ymax):
        return np.ma.filled(self.grid.variables['h'][0][ymin:ymax,xmin:xmax], fill_value=np.nan)

    def read_variable_sea_surface_height_above_mean_sea_level_at_time(self, index_t, xmin, xmax, ymin, ymax):
        return np.ma.filled(self.ncfile.variables["ssh"][index_t][ymin:ymax, xmin:xmax], fill_value=np.nan);

    def read_variable_sea_water_column_thickness_at_time(self, index_t, xmin, xmax, ymin, ymax):
        return np.ma.filled(self.grid.variables["h"][0][ymin:ymax, xmin:xmax] + self.ncfile.variables["ssh"][index_t][ymin:ymax, xmin:xmax], fill_value=np.nan);

    def read_variable_barotropic_sea_water_velocity_at_time(self, index_t, xmin, xmax, ymin, ymax):
        return [np.ma.filled(self.ncfile.variables["u_sea_water_bar_vel"][index_t][ymin:ymax, xmin:xmax], fill_value=np.nan),
                np.ma.filled(self.ncfile.variables["v_sea_water_bar_vel"][index_t][ymin:ymax, xmin:xmax], fill_value=np.nan)]

    #################
    # WAVES
    # Sea Surface
    #################
    def read_variable_sea_surface_wave_significant_height_at_time(self, index_t, xmin, xmax, ymin, ymax):
        return np.ma.filled(self.ncfile.variables["hs"][index_t][ymin:ymax, xmin:xmax], fill_value=np.nan);

    def read_variable_sea_surface_wave_mean_period_at_time(self, index_t, xmin, xmax, ymin, ymax):
        return np.ma.filled(self.ncfile.variables["wave_mean_period"][index_t][ymin:ymax, xmin:xmax],
                            fill_value=np.nan);

    def read_variable_sea_surface_wave_to_direction_at_time(self, index_t, xmin, xmax, ymin, ymax):
        return np.ma.filled(self.ncfile.variables["wave_to_dir"][index_t][ymin:ymax, xmin:xmax], fill_value=np.nan);

    #################
    # METEO
    # At 10 m
    #################
    def read_variable_eastward_wind_10m(self,xmin,xmax,ymin,ymax):
        return np.ma.filled(self.ncfile.variables["u_wind_10m"][:,ymin:ymax,xmin:xmax], fill_value=np.nan)

    def read_variable_northward_wind_10m(self,xmin,xmax,ymin,ymax):
        return np.ma.filled(self.ncfile.variables["v_wind_10m"][:,ymin:ymax,xmin:xmax], fill_value=np.nan)

    def read_variable_wind_10m_at_time(self, index_t,xmin,xmax,ymin,ymax):
        return [np.ma.filled(self.ncfile.variables["u_wind_10m"][index_t][ymin:ymax,xmin:xmax], fill_value=np.nan),
                np.ma.filled(self.ncfile.variables["v_wind_10m"][index_t][ymin:ymax,xmin:xmax], fill_value=np.nan)]



