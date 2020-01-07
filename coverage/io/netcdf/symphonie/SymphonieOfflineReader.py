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
from coverage.io.netcdf.symphonie.SymphonieReader import SymphonieReader
from netCDF4 import Dataset, MFDataset, num2date
import numpy as np
import os

class SymphonieOfflineReader(SymphonieReader):

    def __init__(self,myGrid, myFile):
        SymphonieReader.__init__(self, myGrid, myFile);

    def read_variable_sea_surface_height_above_mean_sea_level_at_time(self,index_t,xmin,xmax,ymin,ymax):
        return np.ma.filled(self.ncfile.variables["ssh_ib"][index_t][ymin:ymax,xmin:xmax],fill_value=np.nan)
        
