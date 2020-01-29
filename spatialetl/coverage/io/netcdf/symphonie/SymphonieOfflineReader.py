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
import logging

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

        
