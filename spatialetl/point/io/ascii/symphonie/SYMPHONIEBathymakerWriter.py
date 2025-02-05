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
 
from __future__ import division, print_function, absolute_import

import pandas

from spatialetl.point.MultiPoint import MultiPoint
from spatialetl.point.io.MultiPointWriter import MultiPointWriter
from spatialetl.utils.VariableDefinition import VariableDefinition
from spatialetl.utils.logger import logging


class SYMPHONIEBathymakerWriter(MultiPointWriter):

    def __init__(self, s,myFile):
        MultiPointWriter.__init__(self, s, myFile);

        if not isinstance(self.points, MultiPoint):
            raise ValueError("This writer supports only a MultiPoint object")

        index = pandas.Index(range(0,self.points.get_nb_points()))
        self.data = pandas.DataFrame(index=index)

        self.data['longitude'] = self.points.read_axis_x()
        self.data['latitude'] = self.points.read_axis_y()

    def close(self):
        self.data.to_csv(self.filename, sep=' ', index=False, columns=["latitude","longitude","bathymetry"], header=False, encoding='utf-8', na_rep="-9999")

    def write_variable_bathymetry(self):
        logging.info('[DefaultMultiPointWriter] Writing variable \''+str(VariableDefinition.LONG_NAME['bathymetry'])+'\'')
        self.data['bathymetry'] = self.points.read_variable_bathymetry()