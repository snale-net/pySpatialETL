##! /usr/bin/env python2.7
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
from spatialetl.coverage.TimeCoverage import TimeCoverage
from netCDF4 import Dataset, num2date
import numpy as np
from datetime import datetime
from time import strftime
import logging
from osgeo import gdal

class SRTMReader(CoverageReader):

    def __init__(self, myFile):
        CoverageReader.__init__(self, myFile);
        self.file = gdal.Open(self.filename)
        self.y = None
        self.x = None

    # Axis
    def read_axis_x(self):

        if self.x is None:
            width = self.file.RasterXSize

            self.x = np.zeros([width])
            for x in range(0,width):
                self.x[x] = self.pixel2coord(x,0)[0]

        return self.x

    def read_axis_y(self):
        if self.y is None:
            height = self.file.RasterYSize

            self.y = np.zeros([height])
            for y in range(0,height):
                self.y[y] = self.pixel2coord(0,y)[1]

        return self.y

    # Scalar
    def read_variable_2D_sea_binary_mask(self):
        width = self.file.RasterXSize
        height = self.file.RasterYSize

        data = np.zeros([height,width])
        data += 1

        return data

    def read_variable_topography(self):

        #print "[ NO DATA VALUE ] = ", band.GetNoDataValue()
        #print "[ MIN ] = ", band.GetMinimum()
        #print "[ MAX ] = ", band.GetMaximum()
        #print "[ SCALE ] = ", band.GetScale()
        #print "[ UNIT TYPE ] = ", band.GetUnitType()

        band = self.file.GetRasterBand(1)
        topo = band.ReadAsArray()
        return [topo]

    def read_variable_bathymetry(self):
        return self.read_variable_topography();

    def pixel2coord(self,y, x):
        # unravel GDAL affine transform parameters
        c, a, b, f, d, e = self.file.GetGeoTransform()
        """Returns global coordinates to pixel center using base-0 raster index"""
        xp = a * y + b * x + a * 0.5 + b * 0.5 + c
        yp = d * y + e * x + d * 0.5 + e * 0.5 + f
        return(xp, yp)
