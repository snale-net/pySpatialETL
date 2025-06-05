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

import numpy as np
from osgeo import gdal

from spatialetl.coverage.io.coverage_reader import CoverageReader


class ImageReader(CoverageReader):

    def __init__(self, myFile):
        CoverageReader.__init__(self, myFile);
        self.file = gdal.Open(self.filename)
        self.y_size = self.file.RasterYSize
        self.x_size = self.file.RasterXSize
        self.y = None
        self.x = None

    def close(self):
        return

    # Axis
    def is_regular_grid(self):
        return True

    def get_x_size(self):
        return self.x_size

    def get_y_size(self):
        return self.y_size

    def read_axis_x(self,xmin,xmax,ymin,ymax):

        if self.x is None:
            width = self.file.RasterXSize

            self.x = np.zeros([width])
            for x in range(0, width):
                self.x[x] = self.pixel2coord(x, 0)[0]

        return self.x[xmin:xmax]

    def read_axis_y(self,xmin,xmax,ymin,ymax):
        if self.y is None:
            height = self.file.RasterYSize

            self.y = np.zeros([height])
            for y in range(0, height):
                self.y[y] = self.pixel2coord(0, y)[1]

        return self.y[ymin:ymax]

    # Scalar
    def read_variable_2D_sea_binary_mask(self,xmin,xmax,ymin,ymax):
        width = self.file.RasterXSize
        height = self.file.RasterYSize

        data = np.zeros([height, width])
        data += 1

        return data[ymin:ymax,xmin:xmax]

    def read_variable_topography(self,xmin,xmax,ymin,ymax):

        # print "[ NO DATA VALUE ] = ", band.GetNoDataValue()
        # print "[ MIN ] = ", band.GetMinimum()
        # print "[ MAX ] = ", band.GetMaximum()
        # print "[ SCALE ] = ", band.GetScale()
        # print "[ UNIT TYPE ] = ", band.GetUnitType()

        band = self.file.GetRasterBand(1)
        topo = band.ReadAsArray()
        return topo[ymin:ymax,xmin:xmax]

    def read_variable_bathymetry(self,xmin,xmax,ymin,ymax):
        return self.read_variable_topography(xmin,xmax,ymin,ymax);

    def pixel2coord(self, y, x):
        # unravel GDAL affine transform parameters
        c, a, b, f, d, e = self.file.GetGeoTransform()
        """Returns global coordinates to pixel center using base-0 raster index"""
        xp = a * y + b * x + a * 0.5 + b * 0.5 + c
        yp = d * y + e * x + d * 0.5 + e * 0.5 + f
        return (xp, yp)



