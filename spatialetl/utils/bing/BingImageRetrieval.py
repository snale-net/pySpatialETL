##! /usr/bin/env python2.7
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

from tempfile import NamedTemporaryFile

from osgeo import gdal

from spatialetl.utils.bing.core.AerialImageRetrieval import AerialImageRetrieval
from spatialetl.utils.logger import logging


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

