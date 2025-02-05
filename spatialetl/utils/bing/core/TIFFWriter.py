#! /usr/bin/env python3.4
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
from osgeo.gdalconst import *

class TIFFWriter ():

    def __init__(self,pattern,myFile):
        # create the output image
        driver = pattern.GetDriver()
        # print driver


        self.tif = driver.Create(myFile,  pattern.RasterXSize, pattern.RasterYSize, 1, GDT_Float32)

        self.tif.SetGeoTransform(pattern.GetGeoTransform())
        self.tif.SetProjection(pattern.GetProjection())

    def write_variable_band_1(self,data):
        outBand = self.tif.GetRasterBand(1)
        # write the data
        outBand.WriteArray(data, 0, 0)

        outBand.FlushCache()
        outBand.SetNoDataValue(0.0)

