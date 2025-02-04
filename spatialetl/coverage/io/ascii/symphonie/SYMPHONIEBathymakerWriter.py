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

import numpy as np

from spatialetl.coverage.io.CoverageWriter import CoverageWriter


class SYMPHONIEBathymakerWriter(CoverageWriter):

    def __init__(self,cov,myFile):
        CoverageWriter.__init__(self,cov,myFile);

        self.file = open(self.filename, "w")

    def write_variable_axis(self):

        lon = self.coverage.read_axis_x()
        lat = self.coverage.read_axis_y()

        file = open(self.filename, "w")
        for i in range(0, len(lon[1])):
            for j in range(0, len(lon)):
                file.write(str(lat[j]) + "\t" + str(lon[i]) + "\t" + "\n")

        file.close()

    def close(self):
        self.file.close()
        
    def write_variable_bathymetry(self): 
        
        lon = self.coverage.read_axis_x()
        lat = self.coverage.read_axis_y()
        data =  self.coverage.read_variable_bathymetry()
        for i in range(0, len(lon)):
            for j in range(0, len(lat)):
                if data[j,i] != np.nan and str(data[j,i]) != "nan" :
                    self.file.write(str(lat[j])+"\t"+str(lon[i])+"\t"+str(data[j,i])+"\n")


