#! /usr/bin/env python2.7
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
import os

class CoverageReader(object):
    
    def __init__(self, myFile):

        if myFile is not None and not os.path.isfile(myFile):
            raise IOError(str(myFile) + " doesn't exists. Abort")

        self.filename = myFile;

    def read_axis_x(self):
        raise NotImplementedError(str(type(self))+" don't have implemented the function 'read_axis_x()'.")

    def read_axis_y(self):
        raise NotImplementedError(str(type(self)) + " don't have implemented the function 'read_axis_y()'.")

    def read_axis_z(self):
        raise NotImplementedError(str(type(self)) + " don't have implemented the function 'read_axis_z()'.")

    def read_axis_t(self, timestamp):
        raise NotImplementedError(str(type(self)) + " don't have implemented the function 'read_axis_t()'.")
