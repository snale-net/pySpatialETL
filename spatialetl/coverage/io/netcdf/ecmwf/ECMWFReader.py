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

import os

import numpy as np
from netCDF4 import Dataset, MFDataset, num2date

from spatialetl.coverage.TimeCoverage import TimeCoverage
from spatialetl.coverage.io.CoverageReader import CoverageReader


class ECMWFReader (CoverageReader):

    def __init__(self, myFile):
        CoverageReader.__init__(self,myFile);

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
        return np.shape(self.ncfile.variables['lon'][:])[0];

    def get_y_size(self):
        return np.shape(self.ncfile.variables['lat'][:])[0];

    def get_t_size(self):
        return np.shape(self.ncfile.variables['time'])[0];

    def read_axis_x(self,xmin,xmax,ymin,ymax):
        return self.ncfile.variables['lon'][xmin:xmax]

    def read_axis_y(self,xmin,xmax,ymin,ymax):
        return self.ncfile.variables['lat'][ymin:ymax]

    def read_axis_t(self,tmin,tmax,timestamp):
        data = self.ncfile.variables['time'][tmin:tmax]
        result = num2date(data, units=self.ncfile.variables['time'].units, calendar="gregorian")

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

    def read_variable_2D_land_binary_mask(self,xmin,xmax,ymin,ymax):
        mask = np.ma.filled(self.ncfile.variables["LSM"][0,ymin:ymax,xmin:xmax], fill_value=np.nan)
        return mask

    def read_variable_2D_sea_binary_mask(self,xmin,xmax,ymin,ymax):
        mask = np.ma.filled(self.ncfile.variables["LSM"][0,ymin:ymax,xmin:xmax], fill_value=np.nan)
        mask += 1.0  # inverse le mask
        mask %= 2  # inverse le mask
        return mask

    def read_variable_3D_sea_binary_mask_at_time(self, index_t,xmin,xmax,ymin,ymax):
        mask = np.ma.filled(self.ncfile.variables["LSM"][index_t,ymin:ymax,xmin:xmax], fill_value=np.nan)
        mask += 1.0  # inverse le mask
        mask %= 2  # inverse le mask
        return mask

    def read_variable_3D_land_binary_mask_at_time(self, index_t,xmin,xmax,ymin,ymax):
        mask = np.ma.filled(self.ncfile.variables["LSM"][index_t,ymin:ymax,xmin:xmax], fill_value=np.nan)
        return mask

    #################
    # METEO
    # 2D
    #################

    def read_variable_rainfall_amount_at_time(self, index_t,xmin,xmax,ymin,ymax):
        return np.ma.filled(self.ncfile.variables["TP"][index_t][ymin:ymax,xmin:xmax], fill_value=np.nan)

    #################
    # METEO
    # Sea surface
    #################

    def read_variable_surface_air_pressure_at_time(self, index_t,xmin,xmax,ymin,ymax):
        sp = np.ma.filled(self.ncfile.variables["SP"][index_t,ymin:ymax,xmin:xmax], fill_value=np.nan)
        sp *= 0.01  # Pa to hPa
        return sp

    def read_variable_sea_surface_air_pressure_at_time(self,index_t,xmin,xmax,ymin,ymax):
        sp = np.ma.filled(self.ncfile.variables["MSL"][index_t,ymin:ymax,xmin:xmax], fill_value=np.nan)
        sp *= 0.01  # Pa to hPa
        return sp

    def read_variable_surface_downward_sensible_heat_flux_at_time(self, index_t,xmin,xmax,ymin,ymax):
        return np.ma.filled(self.ncfile.variables["SSHF"][index_t,ymin:ymax,xmin:xmax], fill_value=np.nan)

    def read_variable_surface_downward_latent_heat_flux_at_time(self, index_t,xmin,xmax,ymin,ymax):
        return np.ma.filled(self.ncfile.variables["SLHF"][index_t,ymin:ymax,xmin:xmax], fill_value=np.nan)

    def read_variable_surface_air_temperature_at_time(self, index_t,xmin,xmax,ymin,ymax):
        return np.ma.filled(self.ncfile.variables["T2M"][index_t,ymin:ymax,xmin:xmax], fill_value=np.nan)

    def read_variable_dew_point_temperature_at_time(self, index_t,xmin,xmax,ymin,ymax):
        return np.ma.filled(self.ncfile.variables["D2M"][index_t,ymin:ymax,xmin:xmax], fill_value=np.nan)

    def read_variable_surface_downwards_solar_radiation_at_time(self, index_t,xmin,xmax,ymin,ymax):
        return np.ma.filled(self.ncfile.variables["SSRD"][index_t,ymin:ymax,xmin:xmax], fill_value=np.nan)

    def read_variable_surface_downwards_thermal_radiation_at_time(self, index_t,xmin,xmax,ymin,ymax):
        return np.ma.filled(self.ncfile.variables["STRD"][index_t,ymin:ymax,xmin:xmax], fill_value=np.nan)

    def read_variable_surface_solar_radiation_at_time(self, index_t,xmin,xmax,ymin,ymax):
        return np.ma.filled(self.ncfile.variables["SSR"][index_t,ymin:ymax,xmin:xmax], fill_value=np.nan)

    def read_variable_surface_thermal_radiation_at_time(self, index_t,xmin,xmax,ymin,ymax):
        return np.ma.filled(self.ncfile.variables["STR"][index_t,ymin:ymax,xmin:xmax], fill_value=np.nan)

    #################
    # METEO
    # At 10 m
    #################

    def read_variable_wind_10m_at_time(self, index_t,xmin,xmax,ymin,ymax):
        return [np.ma.filled(self.ncfile.variables["U10M"][index_t,ymin:ymax,xmin:xmax], fill_value=np.nan), np.ma.filled(self.ncfile.variables["V10M"][index_t,ymin:ymax,xmin:xmax], fill_value=np.nan)]
