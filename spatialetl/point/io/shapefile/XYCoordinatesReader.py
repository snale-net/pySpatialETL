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
#
from __future__ import division, print_function, absolute_import

import numpy as np
from osgeo import ogr

from spatialetl.point.io.MultiPointReader import MultiPointReader


class XYCoordinatesReader(MultiPointReader):
    def __init__(self, myFilename):
        MultiPointReader.__init__(self,myFilename)

        self.shp = ogr.Open(self.filename,0)
        self.lon = []
        self.lat = []
        self.name = []

        if self.shp is None:
            print("pas bon")
        else:
            for layer in self.shp:
                fields = [x.GetName() for x in layer.schema]
                count = 0

                feature = layer.GetNextFeature()
                while feature is not None:

                    flddata = [feature.GetField(feature.GetFieldIndex(x)) for x in fields]
                    g = feature.geometry()
                    attributes = dict(zip(fields, flddata))
                    #attributes["ShpName"] = layer.GetName()
                    if g.GetGeometryType() == 1:  # point
                        self.lon.append(g.GetPoint_2D(0)[0])
                        self.lat.append(g.GetPoint_2D(0)[1])

                        if attributes[fields[0]] is None:
                            self.name.append('Point '+str(count))
                        else:
                            self.name.append(attributes[flddata[0]])
                    # if g.GetGeometryType() == 2:  # linestring
                    # last = g.GetPointCount() - 1
                    # if simplify:
                    #     attributes["Wkb"] = g.ExportToWkb()
                    #     attributes["Wkt"] = g.ExportToWkt()
                    #     attributes["Json"] = g.ExportToJson()
                    #     net.add_edge(g.GetPoint_2D(0), g.GetPoint_2D(last), attributes)
                    # else:
                    #     # separate out each segment as individual edge
                    #     for i in range(last):
                    #         pt1 = g.GetPoint_2D(i)
                    #         pt2 = g.GetPoint_2D(i + 1)
                    #         segment = ogr.Geometry(ogr.wkbLineString)
                    #         segment.AddPoint_2D(pt1[0], pt1[1])
                    #         segment.AddPoint_2D(pt2[0], pt2[1])
                    #         attributes["Wkb"] = segment.ExportToWkb()
                    #         attributes["Wkt"] = segment.ExportToWkt()
                    #         attributes["Json"] = segment.ExportToJson()
                    #         net.add_edge(pt1, pt2, attributes)

                    count = count +1
                    feature = layer.GetNextFeature()

    def read_axis_x(self):
       return self.lon

    def read_axis_y(self):
        return self.lat

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