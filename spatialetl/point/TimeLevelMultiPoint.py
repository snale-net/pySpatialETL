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

from spatialetl.point.LevelMultiPoint import LevelMultiPoint
from spatialetl.point.TimeMultiPoint import TimeMultiPoint


class TimeLevelMultiPoint(LevelMultiPoint, TimeMultiPoint):
    """"""

    def __init__(self,myReader,zbox=None,resolution_z=None,start_time=None,end_time=None,freq=None,time_range=None):
        LevelMultiPoint.__init__(self, myReader,zbox=zbox,resolution_z=resolution_z)
        TimeMultiPoint.__init__(self, myReader,start_time=start_time,end_time=end_time,freq=freq,time_range=time_range)

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

                all_data = self.reader.read_variable_baroclinic_sea_water_velocity_at_time_and_depth(self.map_mpi[self.rank]["src_global_t"].start +indexes_t[t],
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
                    layers[t,z] = self.reader.read_variable_sea_water_temperature_at_time_and_depth(self.map_mpi[self.rank]["src_global_t"].start +indexes_t[t], indexes_z[z])

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
                layers[t, z] = self.reader.read_variable_sea_water_salinity_at_time_and_depth(self.map_mpi[self.rank]["src_global_t"].start +indexes_t[t],
                                                                                                 indexes_z[z])

            results[t] = self.interpolate_vertical(depth, vert_coord, indexes_z, layers[t])

        if len(indexes_t) > 1:
            results = self.interpolate_time(time, indexes_t, results)

        return results


