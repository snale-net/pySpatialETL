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
from spatialetl.coverage.io.netcdf.symphonie.SymphonieReader import SymphonieReader
from spatialetl.utils.VariableDefinition import VariableDefinition
from spatialetl.exception.VariableNameError import VariableNameError
from netCDF4 import Dataset, MFDataset, num2date
import numpy as np
import os
from spatialetl.utils.logger import logging

class SymphonieOfflineReader(SymphonieReader):

    def __init__(self,myGrid, myFile):
        SymphonieReader.__init__(self, myGrid, myFile);

    def read_variable_sea_surface_height_above_mean_sea_level_at_time(self,index_t,xmin,xmax,ymin,ymax):
        try:
            if "ssh_ib" in self.ncfile.variables:
                return np.ma.filled(self.ncfile.variables["ssh_ib"][index_t][ymin:ymax, xmin:xmax], fill_value=np.nan)
        except Exception as ex:
            logging.debug("Error '" + str(ex) + "'")
            raise (VariableNameError("SymphonieOfflineReader", "An error occured : '" + str(ex) + "'", 1000))

        logging.debug("No variables found for '" + str(VariableDefinition.LONG_NAME['sea_surface_height_above_mean_sea_level']) + "'")
        raise (VariableNameError("SymphonieOfflineReader",
                                 "No variables found for '" + str(
                                     VariableDefinition.LONG_NAME['sea_surface_height_above_mean_sea_level']) + "'",
                                 1000))

    def read_variable_barotropic_sea_water_velocity_at_time(self, index_t, xmin, xmax, ymin, ymax):
        raise (VariableNameError("SymphonieOfflineReader",
                                 "No variables found for \'Barotropic Sea Water Velocity\'",
                                 1000))

    def read_variable_baroclinic_sea_water_velocity_at_time_and_depth(self, index_t, index_z, xmin, xmax, ymin, ymax):
        try:
            xmin_overlap, xmax_overlap, ymin_overlap, ymax_overlap, new_xmin, new_xmax, new_ymin, new_ymax = self.compute_overlap_indexes(
                xmin, xmax, ymin, ymax)

            mask_t = self.grid.variables["mask_t"][index_z, ymin_overlap:ymax_overlap, xmin_overlap:xmax_overlap];
            mask_u = self.grid.variables["mask_u"][index_z, ymin_overlap:ymax_overlap, xmin_overlap:xmax_overlap];
            mask_v = self.grid.variables["mask_v"][index_z, ymin_overlap:ymax_overlap, xmin_overlap:xmax_overlap];

            if "u" in self.ncfile.variables:
                data_u = np.ma.filled(
                    self.ncfile.variables["u"][index_t, index_z, ymin_overlap:ymax_overlap,
                    xmin_overlap:xmax_overlap],
                    fill_value=np.nan)
            if "v" in self.ncfile.variables:
                data_v = np.ma.filled(
                    self.ncfile.variables["v"][index_t, index_z, ymin_overlap:ymax_overlap,
                    xmin_overlap:xmax_overlap],
                    fill_value=np.nan)

            u_t, v_t = self.compute_to_tracer(data_u, data_v, mask_t, mask_u, mask_v)

            return [u_t[new_ymin:new_ymax, new_xmin:new_xmax], v_t[new_ymin:new_ymax, new_xmin:new_xmax]]

        except Exception as ex:
            logging.debug("Error '" + str(ex) + "'")
            raise (VariableNameError("SymphonieReader", "An error occured : '" + str(ex) + "'", 1000))

        logging.debug("No variables found for \'Baroclinic Sea Water Velocity\'")
        raise (VariableNameError("SymphonieReader",
                                 "No variables found for \'Baroclinic Sea Water Velocity\'",
                                 1000))


        
