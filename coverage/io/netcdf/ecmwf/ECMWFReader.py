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
from coverage.io.CoverageReader import CoverageReader
from coverage.TimeCoverage import TimeCoverage
from netCDF4 import Dataset,MFDataset, num2date
import numpy as np
from datetime import datetime
from time import strftime
import logging
import os

class ECMWFReader (CoverageReader):

    def __init__(self, myFile):
        CoverageReader.__init__(self,myFile);
        if os.path.isfile(self.filename):
            self.ncfile = Dataset(self.filename, 'r')
        elif os.path.isdir(self.filename):
            self.ncfile = MFDataset(os.path.join(self.filename,"*.nc"), 'r')

    def read_axis_x(self):
        return self.ncfile.variables['lon'][:]

    def read_axis_y(self):
        return self.ncfile.variables['lat'][:]

    def read_axis_t(self, timestamp):
        data = self.ncfile.variables['time'][:]
        result = num2date(data, units=self.ncfile.variables['time'].units, calendar="gregorian")

        if timestamp == 1:
            return [(t - TimeCoverage.TIME_DATUM).total_seconds() \
                    for t in result];
        else:
            return result;

    # Variables

    def read_variable_longitude(self):
        return self.read_axis_x()

    def read_variable_latitude(self):
        return self.read_axis_y()

    def read_variable_time(self):
        return self.read_axis_t(timestamp=0)

    def read_variable_2D_land_binary_mask(self):
        mask = self.ncfile.variables["LSM"][0][:]
        return mask

    def read_variable_2D_sea_binary_mask(self):
        mask = self.ncfile.variables["LSM"][0][:]
        mask += 1.0  # inverse le mask
        mask %= 2  # inverse le mask
        return mask

    def read_variable_3D_sea_binary_mask_at_time(self, index_t):
        mask = self.ncfile.variables["LSM"][index_t][:]
        mask += 1.0  # inverse le mask
        mask %= 2  # inverse le mask
        return mask

    def read_variable_3D_land_binary_mask_at_time(self, index_t):
        mask = self.ncfile.variables["LSM"][index_t][:]
        return mask

    #################
    # METEO
    # 2D
    #################

    def read_variable_rainfall_amount_at_time(self, index_t):
        return self.ncfile.variables["TP"][index_t][:]

    #################
    # METEO
    # Sea surface
    #################

    def read_variable_surface_air_pressure_at_time(self, index_t):
        sp = self.ncfile.variables["SP"][index_t][:]
        sp *= 0.01  # Pa to hPa
        return sp

    def read_variable_sea_surface_air_pressure_at_time(self,index_t):
        sp = self.ncfile.variables["MSL"][index_t][:]
        sp *= 0.01  # Pa to hPa
        return sp

    def read_variable_surface_downward_sensible_heat_flux_at_time(self, index_t):
        return self.ncfile.variables["SSHF"][index_t][:]

    def read_variable_surface_downward_latent_heat_flux_at_time(self, index_t):
        return self.ncfile.variables["SLHF"][index_t][:]

    def read_variable_surface_air_temperature_at_time(self, index_t):
        return self.ncfile.variables["T2M"][index_t][:]

    def read_variable_dew_point_temperature_at_time(self, index_t):
        return self.ncfile.variables["D2M"][index_t][:]

    def read_variable_surface_downwards_solar_radiation_at_time(self, index_t):
        return self.ncfile.variables["SSRD"][index_t][:]

    def read_variable_surface_downwards_thermal_radiation_at_time(self, index_t):
        return self.ncfile.variables["STRD"][index_t][:]

    def read_variable_surface_solar_radiation_at_time(self, index_t):
        return self.ncfile.variables["SSR"][index_t][:]

    def read_variable_surface_thermal_radiation_at_time(self, index_t):
        return self.ncfile.variables["STR"][index_t][:]

    #################
    # METEO
    # At 10 m
    #################

    def read_variable_wind_10m_at_time(self, index_t):
        return [self.ncfile.variables["U10M"][index_t][:], self.ncfile.variables["V10M"][index_t][:]]
