#! /usr/bin/env python3.4
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

import numpy as np
from datetime import datetime
from time import strftime
import logging
from osgeo import gdal

class TIFFReader ():

    def __init__(self, myFile):

        self.file = gdal.Open(myFile)
        self.y = None
        self.x = None

    # Axis
    def read_axis_t(self,timestamp):
        band = self.file.GetRasterBand(1)
        metadata = band.GetMetadata()

        temp = [datetime.fromtimestamp(int(metadata["GRIB_VALID_TIME"].split( )[0]))]

        result = [ datetime.strptime(str(t), '%Y-%m-%d %H:%M:%S') \
                for t in temp];

        return result

    def read_axis_x(self):

        if self.x is None:
            width = self.file.RasterXSize

            self.x = np.zeros([width])
            for x in range(0,width):
                self.x[x] = self.pixel2coord(x,0)[0]
            #print self.x

        return self.x

    def read_axis_y(self):
        if self.y is None:
            height = self.file.RasterYSize

            self.y = np.zeros([height])
            for y in range(0,height):
                self.y[y] = self.pixel2coord(0,y)[1]
            #print self.y

        return self.y

    def getGeoTransform(self):
        return self.file.GetGeoTransform()

    def getProjection(self):
        return self.file.GetProjection()

    # Scalar
    def read_variable_band_1(self):

        #print "[ NO DATA VALUE ] = ", band.GetNoDataValue()
        #print "[ MIN ] = ", band.GetMinimum()
        #print "[ MAX ] = ", band.GetMaximum()
        #print "[ SCALE ] = ", band.GetScale()
        #print "[ UNIT TYPE ] = ", band.GetUnitType()

        band = self.file.GetRasterBand(1)
        return band.ReadAsArray()

    def get_no_data_value_band_1(self):

        band = self.file.GetRasterBand(1)
        return band.GetNoDataValue()

    def pixel2coord(self,y, x):
        # unravel GDAL affine transform parameters
        c, a, b, f, d, e = self.file.GetGeoTransform()
        """Returns global coordinates to pixel center using base-0 raster index"""
        xp = a * y + b * x + a * 0.5 + b * 0.5 + c
        yp = d * y + e * x + d * 0.5 + e * 0.5 + f
        return(xp, yp)
