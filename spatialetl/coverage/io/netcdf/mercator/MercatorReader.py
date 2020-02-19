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
from spatialetl.coverage.io.CoverageReader import CoverageReader
from spatialetl.utils.VariableDefinition import VariableDefinition
from spatialetl.coverage.TimeCoverage import TimeCoverage
from netCDF4 import Dataset, MFDataset, num2date
from spatialetl.exception.VariableNameError import VariableNameError
import numpy as np
import os
from spatialetl.utils.logger import logging

class MercatorReader(CoverageReader):
    
    def __init__(self,myFile):
        CoverageReader.__init__(self,myFile)

        if os.path.isfile(self.filename):
            self.ncfile = Dataset(self.filename, 'r')
        elif os.path.isdir(self.filename):
            self.ncfile = MFDataset(os.path.join(self.filename,"*.nc"), 'r')
        elif self.filename.endswith("*"):
            self.ncfile = MFDataset(self.filename+".nc", 'r')
        else:
            raise ValueError("Unable to decode file "+str(self.filename))

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

    def read_axis_x(self, xmin, xmax, ymin, ymax):
        return self.ncfile.variables['longitude'][xmin:xmax]

    def read_axis_y(self, xmin, xmax, ymin, ymax):
        return self.ncfile.variables['latitude'][ymin:ymax]

    def read_axis_z(self, ):
        lev = self.ncfile.variables["depth"][:]
        # lev = np.ma.filled(self.grid.variables["depth_t"], fill_value=np.nan)
        # lev = np.ma.filled(mx, fill_value=np.nan)
        return lev

    def read_axis_t(self, tmin, tmax, timestamp):
        data = self.ncfile.variables['time'][tmin:tmax]
        result = num2date(data, units=self.ncfile.variables['time'].units, calendar=self.ncfile.variables['time'].calendar)

        if timestamp == 1:
            return [(t - TimeCoverage.TIME_DATUM).total_seconds() \
                    for t in result];
        else:
            return result

    def read_variable_2D_sea_binary_mask(self, xmin, xmax, ymin, ymax):
        try:
            if "zos" in self.ncfile.variables:
                mask = np.ma.filled(self.ncfile.variables["zos"][0,ymin:ymax, xmin:xmax], fill_value=np.nan)
                mask[mask != np.nan]=1
                return mask
        except Exception as ex:
            logging.debug("Error '" + str(ex) + "'")
            raise (VariableNameError("MercatorReader", "An error occured : '" + str(ex) + "'", 1000))

        logging.debug("No variables found for '" + str(VariableDefinition.LONG_NAME['2d_sea_binary_mask']) + "'")
        raise (VariableNameError("MercatorReader",
                                 "No variables found for '" + str(
                                     VariableDefinition.LONG_NAME['2d_sea_binary_mask']) + "'",
                                 1000))

    #################
    # HYDRO
    # Sea Surface
    #################

    def read_variable_sea_surface_height_above_geoid_at_time(self, index_t, xmin, xmax, ymin, ymax):
        try:
            if "zos" in self.ncfile.variables:
                return np.ma.filled(self.ncfile.variables["zos"][index_t,ymin:ymax, xmin:xmax], fill_value=np.nan)
        except Exception as ex:
            logging.debug("Error '" + str(ex) + "'")
            raise (VariableNameError("MercatorReader", "An error occured : '" + str(ex) + "'", 1000))

        logging.debug("No variables found for '" + str(
            VariableDefinition.LONG_NAME['sea_surface_height_above_geoid']) + "'")
        raise (VariableNameError("MercatorReader",
                                 "No variables found for '" + str(
                                     VariableDefinition.LONG_NAME['sea_surface_height_above_geoid']) + "'",
                                 1000))

        #################
        # HYDRO
        # 3D
        #################

    def read_variable_sea_water_temperature_at_time_and_depth(self, index_t, index_z, xmin, xmax, ymin, ymax):
        try:
            if "thetao" in self.ncfile.variables:
                return np.ma.filled(self.ncfile.variables["thetao"][index_t,index_z,ymin:ymax, xmin:xmax],
                                    fill_value=np.nan)
        except Exception as ex:
            logging.debug("Error '" + str(ex) + "'")
            raise (VariableNameError("SymphonieReader", "An error occured : '" + str(ex) + "'", 1000))

        logging.debug("No variables found for '" + str(
            VariableDefinition.LONG_NAME['sea_surface_temperature']) + "'")
        raise (VariableNameError("SymphonieReader",
                                 "No variables found for '" + str(
                                     VariableDefinition.LONG_NAME['sea_surface_temperature']) + "'",
                                 1000))

    def read_variable_sea_water_salinity_at_time_and_depth(self, index_t, index_z, xmin, xmax, ymin, ymax):
        try:
            if "so" in self.ncfile.variables:
                return np.ma.filled(self.ncfile.variables["so"][index_t,index_z,ymin:ymax, xmin:xmax],
                                    fill_value=np.nan)
        except Exception as ex:
            logging.debug("Error '" + str(ex) + "'")
            raise (VariableNameError("SymphonieReader", "An error occured : '" + str(ex) + "'", 1000))

        logging.debug("No variables found for '" + str(
            VariableDefinition.LONG_NAME['sea_surface_salinity']) + "'")
        raise (VariableNameError("SymphonieReader",
                                 "No variables found for '" + str(
                                     VariableDefinition.LONG_NAME['sea_surface_salinity']) + "'",
                                 1000))

    def read_variable_baroclinic_sea_water_velocity_at_time_and_depth(self, index_t, index_z, xmin, xmax, ymin, ymax):
        return [self.ncfile.variables["uo"][index_t,index_z,ymin:ymax, xmin:xmax], self.ncfile.variables["vo"][index_t,index_z,ymin:ymax, xmin:xmax]]
