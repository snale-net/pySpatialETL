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

from spatialetl.coverage.io.CoverageWriter import CoverageWriter


class GMTWriter(CoverageWriter):

    def __init__(self,cov,myFile):
        #raise NotImplementedError("Update to do")
        CoverageWriter.__init__(self,cov,myFile);

        self.file = open(self.filename, "w")

    def write_variable_axis(self):

        lon = self.coverage.read_axis_x()
        lat = self.coverage.read_axis_y()

        file = open(self.filename, "w")
        file.write("#Longitude \t Latitude\n")
        for i in range(0, len(lon[1])):
            for j in range(0, len(lon)):
                # print i,j

                file.write(str(lon[j, i]) + "\t" + str(lat[j, i]) + "\t" + "\n")

        file.close()

    def close(self):
        self.file.close()
        
    def write_variable_bathymetry(self): 
        
        lon = self.coverage.read_axis_x()
        lat = self.coverage.read_axis_y()
        data =  self.coverage.read_variable_bathymetry()
        print(np.shape(data))

        self.file.write("#Longitude \t Latitude \t h (m)\n")
        for i in range(0, len(lon)):
            for j in range(0, len(lat)):
                #print i,j
                
                self.file.write(str(lon[i])+"\t"+str(lat[j])+"\t"+str(data[j,i])+"\n")

        
    def write_variable_current_at_time_and_depth(self,time,z):
        
        lon = self.coverage.read_axis_x()
        lat = self.coverage.read_axis_y()
        cur = self.coverage.read_variable_current_at_time_and_depth(time,z)
        mask = self.coverage.read_variable_2D_mask()
        
        file = open(self.filename, "w")  
        file.write("#Longitude \t Latitude \t u comp (m/s) \t v comp (m/s)\n")
        for i in range(0, self.coverage.get_x_size(),6):
            for j in range(0, self.coverage.get_y_size(),6):
                    if(mask[j,i]==1.):
                        file.write(str(lon[j,i])+"\t"+str(lat[j,i])+"\t"+str(cur[0][j,i])+"\t"+str(cur[1][j,i])+"\n")
                
        file.close()
        
    def write_variable_current_at_time(self,time):
        
        lon = self.coverage.read_axis_x()
        lat = self.coverage.read_axis_y()
        ucur = self.coverage.read_variable_u_current_at_time(time)
        vcur = self.coverage.read_variable_v_current_at_time(time)
        
        file = open(self.filename, "w")  
        file.write("#Longitude \t Latitude \t u comp (m/s) \t v comp (m/s)\n")  
        for i in range(0, self.coverage.get_x_size()):
            for j in range(0, self.coverage.get_y_size()):
                
                file.write(str(lon[j,i])+"\t"+str(lat[j,i])+"\t"+str(ucur[j,i])+"\t"+str(vcur[j,i])+"\n")  
                
        file.close()

    def write_variable_ssh_at_time(self, time):

        lon = self.coverage.read_axis_x()
        lat = self.coverage.read_axis_y()
        ssh = self.coverage.read_variable_ssh_at_time(time)

        file = open(self.filename, "w")
        file.write("#Longitude \t Latitude \t u comp (m/s) \t v comp (m/s)\n")
        for i in range(0, self.coverage.get_x_size()):
            for j in range(0, self.coverage.get_y_size()):
                file.write(
                    str(lon[j, i]) + "\t" + str(lat[j, i]) + "\t" + str(ssh[j, i]) + "\n")

        file.close()

