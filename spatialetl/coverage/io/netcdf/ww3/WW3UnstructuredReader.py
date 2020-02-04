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
from spatialetl.coverage.io.CoverageReader import CoverageReader
from spatialetl.coverage.TimeCoverage import TimeCoverage
from netCDF4 import Dataset, num2date
from spatialetl.utils.VariableDefinition import VariableDefinition
from spatialetl.exception.VariableNameError import VariableNameError
import numpy as np
import cftime
import logging

class WW3UnstructuredReader (CoverageReader):

    def __init__(self,myFile,xsize=None, ysize=None,open_boundaries_mask_to_land=False):
        CoverageReader.__init__(self,myFile);
        self.ncfile = Dataset(self.filename, 'r')

        if xsize is None and ysize is None:
            temp =np.reshape(self.ncfile.variables["MAPSTA"][:], (-1, 1))
            self.x_size =  np.shape(temp)[0]
            self.y_size =  np.shape(temp)[1]
        else:
            self.x_size = xsize
            self.y_size = ysize

        self.open_boundaries_mask_to_land=open_boundaries_mask_to_land

    def is_regular_grid(self):
        return False

    def get_x_size(self):
        return self.x_size

    def get_y_size(self):
        return self.y_size

    def get_t_size(self):
        return np.shape(self.ncfile.variables['time'])[0]
        
    # Axis
    def read_axis_t(self, tmin, tmax, timestamp):
        data = self.ncfile.variables['time'][tmin:tmax]
        temp = num2date(data, units=self.ncfile.variables['time'].units, calendar="julian")

        result = [cftime.datetime(t.year, t.month, t.day, t.hour, t.minute, t.second, t.microsecond) \
                  for t in temp];

        if timestamp == 1:
            return [(t - TimeCoverage.TIME_DATUM).total_seconds() \
                    for t in result];
        else:
            return result;
    
    def read_axis_x(self,xmin,xmax,ymin,ymax):
        return np.reshape(self.ncfile.variables['longitude'][:], (self.y_size, self.x_size))[ymin:ymax,xmin:xmax]
    
    def read_axis_y(self,xmin,xmax,ymin,ymax):
        return np.reshape(self.ncfile.variables['latitude'][:], (self.y_size, self.x_size))[ymin:ymax,xmin:xmax]
    
    # Scalar 
    def read_variable_2D_sea_binary_mask(self,xmin,xmax,ymin,ymax):
        data= np.reshape(self.ncfile.variables["MAPSTA"][:], (self.y_size, self.x_size))[ymin:ymax,xmin:xmax]
        # Boundaries conditions
        if self.open_boundaries_mask_to_land:
            data[data == 2] = 0 # To land
        else:
            data[data == 2] = 1  # To sea
        data = np.ma.filled(data, fill_value=0)
        return data
    
    def read_variable_bathymetry(self,xmin,xmax,ymin,ymax):
        return np.reshape(self.ncfile.variables["dpt"][0][:], (self.y_size, self.x_size))[ymin:ymax,xmin:xmax]

    def read_variable_bathymetry_at_time(self,t,xmin,xmax,ymin,ymax):
        return np.reshape(self.ncfile.variables["dpt"][t][:], (self.y_size, self.x_size))[ymin:ymax,xmin:xmax]
    
    def read_variable_sea_surface_height_above_mean_sea_level_at_time(self,t,xmin,xmax,ymin,ymax):
        return np.reshape(self.ncfile.variables["wlv"][t][:], (self.y_size, self.x_size))[ymin:ymax,xmin:xmax]
    
    def read_variable_sea_surface_wave_significant_height_at_time(self,t,xmin,xmax,ymin,ymax):
        try:
            if "hs" in self.ncfile.variables:
                return np.ma.filled(np.reshape(self.ncfile.variables["hs"][t][:], (self.y_size, self.x_size))[ymin:ymax,xmin:xmax], fill_value=np.nan)
                #print(data)
                #print("*****")
                #mask = self.read_variable_2D_sea_binary_mask(xmin,xmax,ymin,ymax)
                #data = np.ma.filled(np.ma.masked_where((mask != 0) & (mask != 1) ,data), fill_value=np.nan)
                #data = np.ma.filled(np.ma.masked_where((mask != 1), data), fill_value=np.nan)
                #print(data)
        except Exception as ex:
            logging.debug("Error '" + str(ex) + "'")
            raise (VariableNameError("WW3UnstructuredReader", "An error occured : '" + str(ex) + "'", 1000))

        logging.debug("No variables found for '" + str(
            VariableDefinition.LONG_NAME['sea_surface_wave_significant_height']) + "'")
        raise (VariableNameError("WW3UnstructuredReader",
                                 "No variables found for '" + str(
                                     VariableDefinition.LONG_NAME['sea_surface_wave_significant_height']) + "'",
                                 1000))

    def read_variable_sea_surface_wave_breaking_height_at_time(self,t,xmin,xmax,ymin,ymax):
        return np.reshape(self.ncfile.variables["wch"][t][:], (self.y_size, self.x_size))[ymin:ymax,xmin:xmax]
    
    def read_variable_sea_surface_wave_from_direction_at_time(self,t,xmin,xmax,ymin,ymax):
         #attention en ce qui concerne les conventions d'angle en meteorologie. celles ci
        # sont appliquées par ww3. en meteo, un vent venant du nord a une direction de 0°,
        # un vent venant de l'est a une direction de 90°. par consequent il faut corriger
        # cette convention en faisant 270°-angle
        return 270.-np.reshape(self.ncfile.variables["dir"][t][:], (self.y_size, self.x_size))[ymin:ymax,xmin:xmax]

    def read_variable_sea_surface_wave_to_direction_at_time(self, t,xmin,xmax,ymin,ymax):
        return np.reshape(self.ncfile.variables["dir"][t][:], (self.y_size, self.x_size))[ymin:ymax,xmin:xmax]
    
    def read_variable_sea_surface_wave_mean_period_at_time(self,t,xmin,xmax,ymin,ymax):
        return np.reshape(self.ncfile.variables["t01"][t][:], (self.y_size, self.x_size))[ymin:ymax,xmin:xmax]
    
    def read_variable_radiation_pressure_bernouilli_head_at_time(self,t,xmin,xmax,ymin,ymax):
        return np.reshape(self.ncfile.variables["bhd"][t][:], (self.y_size, self.x_size))[ymin:ymax,xmin:xmax]

    def read_variable_sea_surface_wave_energy_flux_to_ocean_at_time(self, t,xmin,xmax,ymin,ymax):
        return np.reshape(self.ncfile.variables["foc"][t][:], (self.y_size, self.x_size))[ymin:ymax,xmin:xmax]

    def read_variable_sea_surface_wave_peak_frequency_at_time(self, t,xmin,xmax,ymin,ymax):
        return np.reshape(self.ncfile.variables["fp"][t][:], (self.y_size, self.x_size))[ymin:ymax,xmin:xmax]
    
    # Vector
    def read_variable_barotropic_sea_water_velocity_at_time(self,t,xmin,xmax,ymin,ymax):
        return [np.reshape(self.ncfile.variables["ucur"][t][:], (self.y_size, self.x_size))[ymin:ymax,xmin:xmax],
                np.reshape(self.ncfile.variables["vcur"][t][:], (self.y_size, self.x_size))[ymin:ymax,xmin:xmax]
               ]

    def read_variable_atmosphere_momentum_flux_to_waves_at_time(self,t,xmin,xmax,ymin,ymax):
        return [np.reshape(self.ncfile.variables["utaw"][t][:], (self.y_size, self.x_size))[ymin:ymax,xmin:xmax],
                np.reshape(self.ncfile.variables["vtaw"][t][:], (self.y_size, self.x_size))[ymin:ymax,xmin:xmax] ]
    
    def read_variable_waves_momentum_flux_to_ocean_at_time(self,t,xmin,xmax,ymin,ymax):
        return [np.reshape(self.ncfile.variables["utwo"][t][:], (self.y_size, self.x_size))[ymin:ymax,xmin:xmax],
                np.reshape(self.ncfile.variables["vtwo"][t][:], (self.y_size, self.x_size))[ymin:ymax,xmin:xmax] ]

    def read_variable_sea_surface_wave_stokes_drift_velocity_at_time(self,t,xmin,xmax,ymin,ymax):
        return [np.reshape(self.ncfile.variables["uuss"][t][:], (self.y_size, self.x_size))[ymin:ymax,xmin:xmax],
                np.reshape(self.ncfile.variables["vuss"][t][:], (self.y_size, self.x_size))[ymin:ymax,xmin:xmax] ]

    def read_variable_wind_10m_at_time(self,t,xmin,xmax,ymin,ymax):
        return [np.reshape(self.ncfile.variables["uwnd"][t][:], (self.y_size, self.x_size))[ymin:ymax,xmin:xmax],
                np.reshape(self.ncfile.variables["vwnd"][t][:], (self.y_size, self.x_size))[ymin:ymax,xmin:xmax] ]
           
    