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
import logging
from osgeo import gdal
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

