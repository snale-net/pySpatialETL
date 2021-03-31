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