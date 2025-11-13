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
import numpy as np
from spatialetl.coverage.io.coverage_reader import CoverageReader


class MemoryReader (CoverageReader):

    def __init__(self,x,y,bathy=None):
        self.x = x
        self.y = y
        self.bathy = bathy

    def is_regular_grid(self):
        return True if len(np.shape(self.x)) == 1 else False

    def get_x_size(self):
        return np.shape(self.x)[0];

    def get_y_size(self):
        return np.shape(self.y)[0];

    def read_axis_x(self,xmin,xmax,ymin=None,ymax=None):
        if self.is_regular_grid():
            return self.x[xmin:xmax]
        else:
            return self.x[ymin:ymax,xmin:xmax]

    def read_axis_y(self,xmin=None,xmax=None,ymin=None,ymax=None):
        if self.is_regular_grid():
            return self.y[ymin:ymax]
        else:
            return self.y[ymin:ymax, xmin:xmax]

    # Variables
    def read_variable_longitude(self,xmin,xmax,ymin,ymax):
        return self.read_axis_x(xmin,xmax,ymin,ymax)

    def read_variable_latitude(self,xmin,xmax,ymin,ymax):
        return self.read_axis_y(xmin,xmax,ymin,ymax)

    def read_variable_2D_land_binary_mask(self,xmin,xmax,ymin,ymax):
        raise NotImplementedError(str(type(self)) + " don't have implemented the function 'get_x_size()'.")

    def read_variable_2D_sea_binary_mask(self,xmin,xmax,ymin,ymax):
        raise NotImplementedError(str(type(self))+" don't have implemented the function 'get_x_size()'.")

    def read_variable_3D_sea_binary_mask_at_time(self, index_t,xmin,xmax,ymin,ymax):
        raise NotImplementedError(str(type(self))+" don't have implemented the function 'get_x_size()'.")

    def read_variable_3D_land_binary_mask_at_time(self, index_t,xmin,xmax,ymin,ymax):
        raise NotImplementedError(str(type(self)) + " don't have implemented the function 'get_x_size()'.")

    #################
    # HYDRO
    # 2D
    #################

    def read_variable_bathymetry(self, xmin, xmax, ymin, ymax):
       return self.bathy[ymin:ymax,xmin:xmax]

    #################
    # METEO
    # 2D
    #################

    def read_variable_rainfall_amount_at_time(self, index_t,xmin,xmax,ymin,ymax):
        raise NotImplementedError(str(type(self)) + " don't have implemented the function 'get_x_size()'.")

    #################
    # METEO
    # Sea surface
    #################

    def read_variable_surface_air_pressure_at_time(self, index_t,xmin,xmax,ymin,ymax):
        raise NotImplementedError(str(type(self))+" don't have implemented the function 'get_x_size()'.")

    def read_variable_sea_surface_air_pressure_at_time(self,index_t,xmin,xmax,ymin,ymax):
        raise NotImplementedError(str(type(self))+" don't have implemented the function 'get_x_size()'.")

    def read_variable_surface_downward_sensible_heat_flux_at_time(self, index_t,xmin,xmax,ymin,ymax):
        raise NotImplementedError(str(type(self)) + " don't have implemented the function 'get_x_size()'.")

    def read_variable_surface_downward_latent_heat_flux_at_time(self, index_t,xmin,xmax,ymin,ymax):
        raise NotImplementedError(str(type(self)) + " don't have implemented the function 'get_x_size()'.")

    def read_variable_surface_air_temperature_at_time(self, index_t,xmin,xmax,ymin,ymax):
        raise NotImplementedError(str(type(self))+" don't have implemented the function 'get_x_size()'.")

    def read_variable_dew_point_temperature_at_time(self, index_t,xmin,xmax,ymin,ymax):
        raise NotImplementedError(str(type(self)) + " don't have implemented the function 'get_x_size()'.")

    def read_variable_surface_downward_solar_radiation_at_time(self, index_t, xmin, xmax, ymin, ymax):
        raise NotImplementedError(str(type(self))+" don't have implemented the function 'get_x_size()'.")

    def read_variable_surface_downward_thermal_radiation_at_time(self, index_t, xmin, xmax, ymin, ymax):
        raise NotImplementedError(str(type(self)) + " don't have implemented the function 'get_x_size()'.")

    def read_variable_surface_solar_radiation_at_time(self, index_t,xmin,xmax,ymin,ymax):
        raise NotImplementedError(str(type(self))+" don't have implemented the function 'get_x_size()'.")

    def read_variable_surface_thermal_radiation_at_time(self, index_t,xmin,xmax,ymin,ymax):
        raise NotImplementedError(str(type(self)) + " don't have implemented the function 'get_x_size()'.")

    #################
    # METEO
    # At 10 m
    #################

    def read_variable_wind_10m_at_time(self, index_t,xmin,xmax,ymin,ymax):
        raise NotImplementedError(str(type(self)) + " don't have implemented the function 'get_x_size()'.")
