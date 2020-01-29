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

from spatialetl.point.LevelMultiPoint import LevelMultiPoint
from spatialetl.point.TimeMultiPoint import TimeMultiPoint
from spatialetl.point.MultiPoint import MultiPoint

import numpy as np
from datetime import datetime,timedelta
import logging


class TimeLevelMultiPoint(LevelMultiPoint, TimeMultiPoint):
    """"""

    def __init__(self,myReader,depth=None,time_range=None):
        LevelMultiPoint.__init__(self, myReader,depth)
        TimeMultiPoint.__init__(self, myReader,time_range=time_range)

    # Scalar
    def read_variable_baroclinic_sea_water_velocity_at_time_and_depth(self,time,depth):

        indexes_t = self.find_time_index(time);
        tmp = self.find_level_index(depth);
        vert_coord = tmp[0]
        indexes_z = tmp[1]

        layers = np.zeros([np.shape(indexes_t)[0],2, np.shape(indexes_z)[0],self.get_nb_points()])
        layers[::] = np.NAN

        results = np.zeros([np.shape(indexes_t)[0],2,self.get_nb_points()])
        results[:] = np.NAN

        for t in range(0, len(indexes_t)):

            for z in range(0, len(indexes_z)):

                all_data = self.reader.read_variable_baroclinic_sea_water_velocity_at_time_and_depth(indexes_t[t],
                                                                                                 indexes_z[z])

                # Comp U
                layers[t,0, z] = all_data[0]
                # Comp V
                layers[t,1, z] = all_data[1]

            # Comp U
            results[t][0] = self.interpolate_vertical(depth, vert_coord, indexes_z, layers[t,0])
            # Comp V
            results[t][1] = self.interpolate_vertical(depth, vert_coord, indexes_z, layers[t,1])

        if len(indexes_t) == 1:
            return [results[0:,0,:],results[0:,1,:]]
        else:
            return [self.interpolate_time(time, indexes_t, results[:,0,:]),self.interpolate_time(time, indexes_t, results[:,1,:])]


    def read_variable_sea_water_temperature_at_time_and_depth(self,time,depth):

        indexes_t = self.find_time_index(time);
        tmp = self.find_level_index(depth);
        vert_coord = tmp[0]
        indexes_z = tmp[1]

        layers = np.zeros([np.shape(indexes_t)[0],np.shape(indexes_z)[0], self.get_nb_points()])
        layers[::] = np.NAN

        results = np.zeros([np.shape(indexes_t)[0],self.get_nb_points()])
        results[:] = np.NAN

        for t in range(0, len(indexes_t)):
            for z in range(0, len(indexes_z)):
                    layers[t,z] = self.reader.read_variable_sea_water_temperature_at_time_and_depth(indexes_t[t], indexes_z[z])

            results[t] = self.interpolate_vertical(depth,vert_coord,indexes_z,layers[t])

        if len(indexes_t) > 1:
            results = self.interpolate_time(time, indexes_t, results)

        return results


    def read_variable_sea_water_salinity_at_time_and_depth(self, time, depth):

        indexes_t = self.find_time_index(time);
        tmp = self.find_level_index(depth);
        vert_coord = tmp[0]
        indexes_z = tmp[1]

        layers = np.zeros([np.shape(indexes_t)[0], np.shape(indexes_z)[0], self.get_nb_points()])
        layers[::] = np.NAN

        results = np.zeros([np.shape(indexes_t)[0], self.get_nb_points()])
        results[:] = np.NAN

        for t in range(0, len(indexes_t)):
            for z in range(0, len(indexes_z)):
                layers[t, z] = self.reader.read_variable_sea_water_salinity_at_time_and_depth(indexes_t[t],
                                                                                                 indexes_z[z])

            results[t] = self.interpolate_vertical(depth, vert_coord, indexes_z, layers[t])

        if len(indexes_t) > 1:
            results = self.interpolate_time(time, indexes_t, results)

        return results


