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
from spatialetl.point.io.MultiPointReader import MultiPointReader
import os.path
import pandas
import logging
import math
import numpy as np


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