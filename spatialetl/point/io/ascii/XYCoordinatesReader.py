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

import numpy as np
import pandas

from spatialetl.point.io.MultiPointReader import MultiPointReader


class XYCoordinatesReader(MultiPointReader):
    def __init__(self, myFilename):
        MultiPointReader.__init__(self,myFilename)

        self.data = pandas.read_csv(self.filename, usecols=[0 , 1, 2],
                               names=['longitude','latitude','name'], sep='\t')

    def read_axis_x(self):
        return self.data['longitude'].values

    def read_axis_y(self):
        return self.data['latitude'].values

    def get_point_names(self):
        return self.data['name'].values

    def get_coordinates(self):
        x = self.read_axis_x()
        y = self.read_axis_y()
        nbPoints = np.shape(x)[0]

        data = np.zeros([nbPoints,2])

        for i in range(0,nbPoints):
            data[i] = [x[i],y[i]]

        return data