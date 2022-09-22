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

from netCDF4 import Dataset
import numpy as np

from spatialetl.coverage.io.CoverageReader import CoverageReader


class GMTReader(CoverageReader):

    def __init__(self, myFile):
        CoverageReader.__init__(self,myFile);
        self.ncfile = Dataset(self.filename, 'r')

    def close(self):
        self.ncfile.close()

    def is_regular_grid(self):
        return True

    def get_x_size(self):
        return np.shape(self.ncfile.variables['x'][:])[0];

    def get_y_size(self):
        return np.shape(self.ncfile.variables['y'][:])[0];

    def read_axis_x(self, xmin, xmax, ymin, ymax):
        return self.ncfile.variables['x'][xmin:xmax]

    def read_axis_y(self, xmin, xmax, ymin, ymax):
        return self.ncfile.variables['y'][ymin:ymax]


    # Variables
    def read_variable_longitude(self,xmin, xmax, ymin, ymax):
        return self.read_axis_x(xmin, xmax, ymin, ymax)

    def read_variable_latitude(self,xmin, xmax, ymin, ymax):
        return self.read_axis_y(xmin, xmax, ymin, ymax)

    def read_variable_2D_sea_binary_mask(self,xmin, xmax, ymin, ymax):
        mask = self.ncfile.variables["z"][ymin:ymax,xmin:xmax]
        return mask

    def read_variable_bathymetry(self,xmin, xmax, ymin, ymax):
        return self.ncfile.variables["z"][ymin:ymax,xmin:xmax]
