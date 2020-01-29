#! #! /usr/bin/env python2.7
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
from netCDF4 import Dataset, num2date
import numpy as np
from spatialetl.coverage.io.netcdf.symphonie.SymphonieReader import SymphonieReader

class SymphonieBathycoteReader(SymphonieReader):
    """
    """
    def __init__(self,myFile):
        SymphonieReader.__init__(self,myFile,myFile);
        
    def read_variable_2D_sea_binary_mask(self):
        return self.grid.variables["mask_t"][:]

    def read_variable_mesh_size(self):
        data= self.grid.variables["mesh_size"][:]
        data[data < 0] = np.nan
        return data
