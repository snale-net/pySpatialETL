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
                self.file.write(str(lat[j])+"\t"+str(lon[i])+"\t"+str(data[j,i])+"\n")


