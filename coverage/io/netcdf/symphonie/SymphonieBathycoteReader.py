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
from __future__ import division, print_function, absolute_import
from coverage.io.CoverageReader import CoverageReader
from netCDF4 import Dataset, num2date
import numpy as np

class SymphonieBathycoteReader(CoverageReader):
    """
La classe SymphonieReader permet de lire les données du format Symphonie

@param  myGrid : lien vers le fichier de grille (que l'on trouve dans le RDIR/tmp/grid.nc)
@param myFile : lien vers le fichier de données (que l'on trouve dans GRAPHIQUES)
"""
    def __init__(self,myFile):
        CoverageReader.__init__(self,myFile);
        self.ncfile = Dataset(self.filename, 'r')
        
    # Axis
    def read_axis_t(self,timestamp=0):
        data = self.ncfile.variables['time'][:]         
        result = num2date(data, units = self.ncfile.variables['time'].units.replace('from','since').replace('mar','03').replace('feb','02').replace('jun','06'), calendar = self.ncfile.variables['time'].calendar)
        
        if timestamp ==1:           
            return [ (t - TimeCoverage.TIME_DATUM).total_seconds() \
                for t in result];
        else:            
            return result
    
    def read_axis_x(self):       
        return self.ncfile.variables['longitude_t'][:]
    
    def read_axis_y(self):      
        return self.ncfile.variables['latitude_t'][:]
    
    def read_axis_z(self):
        lev = self.ncfile.variables['depth_t'][::]
        lev[::] *= -1.0 # inverse la profondeur
        return lev
        
    # Data    
    def read_variable_2D_mask(self):
        return self.ncfile.variables["mask_t"][0][:]

    def read_variable_3D_mask(self):
        return self.ncfile.variables["mask_t"][::]
    
    def read_variable_mesh_size(self): 
        data= self.ncfile.variables["mesh_size"][:]
        data[data < 0] = np.nan
        return data
    
    def read_variable_bathymetry(self): 
        return self.ncfile.variables["hm_w"][:]
