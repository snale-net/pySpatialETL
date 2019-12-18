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
from coverage.io.CoverageReader import CoverageReader
from coverage.TimeCoverage import TimeCoverage
from netCDF4 import Dataset, MFDataset, num2date
import numpy as np
from datetime import datetime
from time import strftime
import logging
import os

class WW3Reader (CoverageReader):

    def __init__(self, myFile):
        CoverageReader.__init__(self,myFile);

        if os.path.isfile(self.filename):
            self.ncfile = Dataset(self.filename, 'r')
        elif os.path.isdir(self.filename):
            self.ncfile = MFDataset(os.path.join(self.filename, "*.nc"), 'r')
        
    # Axis
    def read_axis_t(self,timestamp):
        data = self.ncfile.variables['time'][:]        
        temp = num2date(data, units = self.ncfile.variables['time'].units, calendar = "julian")
        
        result = [ datetime.strptime(str(t), '%Y-%m-%d %H:%M:%S') \
                for t in temp];
      
        if timestamp ==1:           
            return [ (t - TimeCoverage.TIME_DATUM).total_seconds() \
                for t in result];
        else:            
            return result;
    
    def read_axis_x(self):        
        return self.ncfile.variables['longitude'][:]
    
    def read_axis_y(self):        
        return self.ncfile.variables['latitude'][:]
    
    # Scalar 
    def read_variable_2D_sea_binary_mask(self):
        mask = self.ncfile.variables["MAPSTA"][:]
        #mask += 1.0 # inverse le mask
        #mask %= 2 # inverse le mask
        return mask

    def read_variable_bathymetry(self):
        return self.ncfile.variables["dpt"][0][:]

    def read_variable_bathymetry_at_time(self,t):
        return self.ncfile.variables["dpt"][t][:]
    
    def read_variable_sea_surface_height_above_mean_sea_level_at_time(self,t):
        return self.ncfile.variables["wlv"][t][:]
    
    def read_variable_sea_surface_wave_significant_height_at_time(self,t):
        return self.ncfile.variables["hs"][t][:]

    def read_variable_sea_surface_wave_breaking_height_at_time(self,t):
        return self.ncfile.variables["wch"][t][:]
    
    def read_variable_sea_surface_wave_from_direction_at_time(self,t):
        #attention en ce qui concerne les conventions d'angle en meteorologie. celles ci
        # sont appliquées par ww3. en meteo, un vent venant du nord a une direction de 0°,
        # un vent venant de l'est a une direction de 90°. par consequent il faut corriger
        # cette convention si on veut la direction vers en faisant 270°-angle
        return self.ncfile.variables["dir"][t][:]

    def read_variable_sea_surface_wave_to_direction_at_time(self, t):
        #attention en ce qui concerne les conventions d'angle en meteorologie. celles ci
        # sont appliquées par ww3. en meteo, un vent venant du nord a une direction de 0°,
        # un vent venant de l'est a une direction de 90°. par consequent il faut corriger
        # cette convention si on veut la direction vers en faisant 270°-angle
        return 270.-self.ncfile.variables["dir"][t][:]
    
    def read_variable_sea_surface_wave_mean_period_at_time(self,t):
        return self.ncfile.variables["t01"][t][:]

    def read_variable_sea_surface_wave_peak_period_at_time(self,t):
        return self.ncfile.variables["tp"][t][:]

    def read_variable_sea_surface_wave_energy_dissipation_at_ground_level_at_time(self,t):
        return self.ncfile.variables["fbb"][t][:]
    
    def read_variable_radiation_pressure_bernouilli_head_at_time(self,t):
        return self.ncfile.variables["bhd"][t][:]

    def read_variable_sea_surface_wave_energy_flux_to_ocean_at_time(self,t):
        return self.ncfile.variables["foc"][t][:]

    def read_variable_sea_surface_wave_peak_frequency_at_time(self,t):
        return self.ncfile.variables["fp"][t][:]
    
    # Vector
    def read_variable_barotropic_sea_water_velocity_at_time(self,t):
        return [self.ncfile.variables["ucur"][t][:], self.ncfile.variables["vcur"][t][:]]

    def read_variable_atmosphere_momentum_flux_to_waves_at_time(self,t):
        return [self.ncfile.variables["utaw"][t][:],self.ncfile.variables["vtaw"][t][:]]

    def read_variable_waves_momentum_flux_to_ocean_at_time(self,t):
        return [self.ncfile.variables["utwo"][t][:],self.ncfile.variables["vtwo"][t][:]]

    def read_variable_sea_surface_wave_stokes_drift_velocity_at_time(self,t):
        return [self.ncfile.variables["uuss"][t][:],self.ncfile.variables["vuss"][t][:]]

    def read_variable_wind_10m_at_time(self,t):
        return [self.ncfile.variables["uwnd"][t][:],self.ncfile.variables["vwnd"][t][:]]
           
    