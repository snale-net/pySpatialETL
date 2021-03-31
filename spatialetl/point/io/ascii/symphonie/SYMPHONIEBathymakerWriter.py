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