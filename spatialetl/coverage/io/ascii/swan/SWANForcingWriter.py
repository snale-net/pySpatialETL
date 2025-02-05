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

from spatialetl.coverage.io.CoverageWriter import CoverageWriter
from spatialetl.utils.logger import logging


class SWANForcingWriter(CoverageWriter):

    def __init__(self,cov,myFile):
        CoverageWriter.__init__(self,cov,myFile);

        self.file = open(self.filename, "w")

    def close(self):
        self.file.close()

    def write_variable_wind_10m(self):

        if self.coverage.rank == 0:
            logging.info('[SWANForcingWriter] Gathering data from all processors')
            global_data = np.empty(
                [self.coverage.get_t_size(type="target_global"),
                 2,
                 self.coverage.get_y_size(type="target_global"),
                 self.coverage.get_x_size(type="target_global")])
            global_data[:] = np.nan

        for time_index in range(0,self.coverage.get_t_size()):

            local_data = self.coverage.read_variable_wind_10m_at_time(time_index)

            if self.coverage.rank != 0:
                self.coverage.comm.Send(np.ascontiguousarray(local_data), dest=0)
            else:
                # Pour le proc n°1
                global_data[self.coverage.map_mpi[self.coverage.rank]["dst_global_t"].start + time_index,
                            0,
                            self.coverage.map_mpi[self.coverage.rank]["dst_global_y"],
                            self.coverage.map_mpi[self.coverage.rank]["dst_global_x"]] = local_data[0]

                global_data[self.coverage.map_mpi[self.coverage.rank]["dst_global_t"].start + time_index,
                            1,
                            self.coverage.map_mpi[self.coverage.rank]["dst_global_y"],
                            self.coverage.map_mpi[self.coverage.rank]["dst_global_x"]] = local_data[1]

                # Pour les autres
                for source in range(1,self.coverage.size):
                    recvbuf = np.empty([2,self.coverage.map_mpi[source]["dst_local_y_size"],self.coverage.map_mpi[source]["dst_local_x_size"]])
                    self.coverage.comm.Recv(recvbuf,source=source)

                    global_data[self.coverage.map_mpi[source]["dst_global_t"].start+time_index,
                                0,
                         self.coverage.map_mpi[source]["dst_global_y"],
                         self.coverage.map_mpi[source]["dst_global_x"]] = recvbuf[0]

                    global_data[self.coverage.map_mpi[source]["dst_global_t"].start + time_index,
                                1,
                                self.coverage.map_mpi[source]["dst_global_y"],
                                self.coverage.map_mpi[source]["dst_global_x"]] = recvbuf[1]

        self.coverage.comm.barrier()

        if self.coverage.rank == 0:

            logging.info('[SWANForcingWriter] Writing variable \'Wind 10m\'')

            file = open(self.filename, "w")

            for time_index in range(0,self.coverage.get_t_size(type="target_global")):
                time = self.coverage.read_axis_t(type="target_global")[time_index]

                logging.debug('[SWANForcingWriter] Writing variable \'Wind 10m\' at time \'' + str(time) + '\'')

                file.write(time.strftime("%Y%m%d.%H%M%S")+"\n")

                for vector in range(0,2):
                    if vector == 0:
                        file.write("u-component\n")
                    else:
                        file.write("v-component\n")
                    for i in range(0, self.coverage.get_x_size()):
                        for j in range(0, self.coverage.get_y_size()):
                            file.write(str(global_data[time_index][vector][j,i])+"\n")

            file.close()




