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
import logging
from tempfile import NamedTemporaryFile
from osgeo import gdal
from spatialetl.utils.bing import AerialImageRetrieval

class BingImageRetrieval():

    def __init__(self,myFile,xmin,xmax,ymin,ymax):
        self.filename = myFile
        jpegFile = NamedTemporaryFile(suffix=".jpeg").name
        f = open(myFile,"w")
        f.close()

        imgretrieval = AerialImageRetrieval(ymin,xmin,ymax,xmax,jpegFile)

        zoom = imgretrieval.max_resolution_imagery_retrieval()
        if zoom != -999:
            gdal.Translate(self.filename,jpegFile,format='GTiff',
                           outputSRS='EPSG:4326', outputBounds=[xmin,ymax,xmax,ymin])

            logging.info("[BingImageRetrieval] Satellite image was successfully retrieved with the maximum resolution available (zoom level {})".format(
                zoom))
        else:
            logging.error("[BingImageRetrieval] Satellite image cannot be retrieved. Possible reason: expected tile image does not exist.")
            raise ValueError("[BingImageRetrieval] Satellite image cannot be retrieved. Possible reason: expected tile image does not exist.")

