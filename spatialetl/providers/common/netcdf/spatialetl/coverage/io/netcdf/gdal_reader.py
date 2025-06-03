#! /usr/bin/env python2.7
# -*- coding: utf-8 -*-
# MIT License
# Copyright (c) 2024 [SNALE - French SAS Company - RCS 951 724 616]
#
# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:
#
# The above copyright notice and this permission notice shall be included in all
# copies or substantial portions of the Software.
#
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
# SOFTWARE.
from __future__ import division, print_function, absolute_import

from netCDF4 import Dataset

from spatialetl.coverage.io.coverage_reader import CoverageReader


class GDALReader (CoverageReader):

    def __init__(self, myFile):
        CoverageReader.__init__(self,myFile);
        self.ncfile = Dataset(self.filename, 'r')

    def close(self):
        self.ncfile.close()

    def is_regular_grid(self):
        return True

    def read_axis_x(self,xmin,xmax,ymin,ymax):
        return self.ncfile.variables['lon'][xmin:xmax]

    def read_axis_y(self,xmin,xmax,ymin,ymax):
        return self.ncfile.variables['lat'][ymin:ymax]

    # Variables
    def read_variable_longitude(self,xmin,xmax,ymin,ymax):
        return self.read_axis_x(xmin,xmax,ymin,ymax)

    def read_variable_latitude(self,xmin,xmax,ymin,ymax):
        return self.read_axis_y(xmin,xmax,ymin,ymax)

    def read_variable_2D_sea_binary_mask(self,xmin,xmax,ymin,ymax):
        mask = self.ncfile.variables["Band1"][:]
        return mask

    def read_variable_bathymetry(self,xmin,xmax,ymin,ymax):
        return self.ncfile.variables["Band1"][ymin:ymax,xmin:xmax]
