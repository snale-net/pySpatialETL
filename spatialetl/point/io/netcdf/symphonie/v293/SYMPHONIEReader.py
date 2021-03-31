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

from spatialetl.coverage.io.netcdf.symphonie.v293.SYMPHONIEReader import SYMPHONIEReader as CovReader
from spatialetl.point.io.netcdf.symphonie.AbstractSYMPHONIEReader import AbstractSYMPHONIEReader


class SYMPHONIEReader(AbstractSYMPHONIEReader):

    def __init__(self,myGrid,myFile,xy,names=None):
        AbstractSYMPHONIEReader.__init__(self, myFile, xy, names);
        self.reader = CovReader(myGrid, self.filename)
        self.find_points_coordinates(xy)



