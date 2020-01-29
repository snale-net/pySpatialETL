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
from spatialetl.coverage.io.CoverageWriter import CoverageWriter
import numpy as np

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

