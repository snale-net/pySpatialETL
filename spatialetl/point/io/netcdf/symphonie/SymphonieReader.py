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
from spatialetl.point.io.MultiPointReader import MultiPointReader
from spatialetl.point.TimeMultiPoint import TimeMultiPoint
from spatialetl.coverage.TimeLevelCoverage import TimeLevelCoverage
from spatialetl.coverage.io.netcdf.symphonie.SymphonieReader import SymphonieReader as CovReader
from netCDF4 import num2date
from scipy.io import loadmat
import numpy as np
from datetime import datetime,timedelta,timezone
import pytz
import logging


class SymphonieReader(MultiPointReader):

    def __init__(self,myGrid,myFile,xy,names=None):
        MultiPointReader.__init__(self, myFile);
        self.coverage = TimeLevelCoverage(CovReader(myGrid,self.filename))
        self.nbPoints = np.shape(xy)[0]
        self.xy_coords = np.zeros([self.nbPoints,2],dtype=np.int32)
        self.xy_values = np.zeros([self.nbPoints, 2])
        self.meta_data = ""

        if names is None:
            self.names = np.empty([self.nbPoints], dtype=object)
            for count in range(0, self.nbPoints):
                self.names[count] = "Point-" + str(count)
        else:
            self.names = names

        for i in range(0, np.shape(self.xy_coords)[0]):
            nearestPoint = self.coverage.find_point_index(xy[i][0], xy[i][1])
            logging.info("Nearest point : " + str(nearestPoint[2]) + " / " + str(nearestPoint[3]) + " at " + str(round(nearestPoint[4], 4)) + " km")
            self.meta_data = self.meta_data + "\n# " + str(self.names[i]) + " : nearest point in SYMPHONIE file is " + str(
                round(nearestPoint[4], 4)) + " km from the target point"
            logging.info("Nearest point (i,j) : " + str(nearestPoint[0]) + " / " + str(nearestPoint[1]))
            self.xy_coords[i] = [nearestPoint[0],nearestPoint[1]]
            self.xy_values[i] = [nearestPoint[2], nearestPoint[3]]

        self.names = names

    # Axis
    def read_axis_x(self):
        return self.xy_values[:,0]

    def read_axis_y(self):
        return self.xy_values[:,1]

    def read_axis_z(self):
        data = np.zeros([self.nbPoints,np.shape(self.coverage.read_axis_z())[0]])
        data[:] = np.nan

        all_data = self.coverage.read_axis_z()

        for index_x in range(0, self.nbPoints):
            data[index_x] = all_data[:,self.xy_coords[index_x][1], self.xy_coords[index_x][0]]

        return data

    def read_axis_t(self,timestamp=0):
        return self.coverage.read_axis_t(timestamp)

    def read_metadata(self):
        m = {}
        m["data_source"] = "SYMPHONIE file"
        m["meta_data"] = self.meta_data

        return m

    def read_variable_longitude(self):
        return self.read_axis_x()

    def read_variable_latitude(self):
        return self.read_axis_y()

    def read_variable_depth(self):
        return self.read_axis_z()

    def read_variable_time(self):
        return self.read_axis_t(timestamp=0)

    def read_variable_point_names(self):
        return self.names

    #################
    # HYDRO
    # Sea Surface
    #################
    def read_variable_sea_surface_height_above_mean_sea_level_at_time(self,index_t):

        data = np.zeros([self.nbPoints])
        data[:] = np.nan

        all_data = self.coverage.reader.read_variable_sea_surface_height_above_mean_sea_level_at_time(index_t)

        for index_x in range(0, self.nbPoints):
            data[index_x] = all_data[self.xy_coords[index_x][1], self.xy_coords[index_x][0]]

        return data

    def read_variable_sea_surface_temperature_at_time(self,index_t):
        data = np.zeros([self.nbPoints])
        data[:] = np.nan

        all_data = self.coverage.reader.read_variable_sea_surface_temperature_at_time(index_t)

        for index_x in range(0, self.nbPoints):
            data[index_x] = all_data[self.xy_coords[index_x][1], self.xy_coords[index_x][0]]

        return data

    def read_variable_sea_surface_salinity_at_time(self,index_t):
        data = np.zeros([self.nbPoints])
        data[:] = np.nan

        all_data = self.coverage.reader.read_variable_sea_surface_salinity_at_time(index_t)

        for index_x in range(0, self.nbPoints):
            data[index_x] = all_data[self.xy_coords[index_x][1], self.xy_coords[index_x][0]]

        return data

    def read_variable_sea_water_velocity_at_sea_water_surface_at_time(self,index_t):

        data = np.zeros([2, self.nbPoints])
        data[:] = np.nan

        all_data = self.coverage.reader.read_variable_sea_water_velocity_at_sea_water_surface_at_time(index_t)

        for index_x in range(0, self.nbPoints):
            data[0][index_x] = all_data[0][self.xy_coords[index_x][1], self.xy_coords[index_x][0]]
            data[1][index_x] = all_data[1][self.xy_coords[index_x][1], self.xy_coords[index_x][0]]

        return data

    #################
    # HYDRO
    # Ground level
    #################
    def read_variable_sea_water_temperature_at_ground_level_at_time(self,index_t):
        data = np.zeros([self.nbPoints])
        data[:] = np.nan

        all_data = self.coverage.reader.read_variable_sea_water_temperature_at_ground_level_at_time(index_t)

        for index_x in range(0, self.nbPoints):
            data[index_x] = all_data[self.xy_coords[index_x][1], self.xy_coords[index_x][0]]

        return data
    def read_variable_sea_water_salinity_at_ground_level_at_time(self,index_t):
        data = np.zeros([self.nbPoints])
        data[:] = np.nan

        all_data = self.coverage.reader.read_variable_sea_water_salinity_at_ground_level_at_time(index_t)

        for index_x in range(0, self.nbPoints):
            data[index_x] = all_data[self.xy_coords[index_x][1], self.xy_coords[index_x][0]]

        return data

    def read_variable_sea_water_velocity_at_ground_level_at_time(self,index_t):

        data = np.zeros([2, self.nbPoints])
        data[:] = np.nan

        all_data = self.coverage.reader.read_variable_sea_water_velocity_at_ground_level_at_time(index_t)

        for index_x in range(0, self.nbPoints):
            data[0][index_x] = all_data[0][self.xy_coords[index_x][1], self.xy_coords[index_x][0]]
            data[1][index_x] = all_data[1][self.xy_coords[index_x][1], self.xy_coords[index_x][0]]

        return data

    #################
    # HYDRO
    # 2D
    #################
    def read_variable_bathymetry(self):
        data = np.zeros([self.nbPoints])
        data[:] = np.nan

        all_data = self.coverage.reader.read_variable_bathymetry()

        for index_x in range(0, self.nbPoints):
            data[index_x] = all_data[self.xy_coords[index_x][1], self.xy_coords[index_x][0]]

        return data

    def read_variable_barotropic_sea_water_velocity_at_time(self,index_t):
        data = np.zeros([2, self.nbPoints])
        data[:] = np.nan

        all_data = self.coverage.reader.read_variable_barotropic_sea_water_velocity_at_time(index_t)

        for index_x in range(0, self.nbPoints):
            data[0][index_x] = all_data[0][self.xy_coords[index_x][1], self.xy_coords[index_x][0]]
            data[1][index_x] = all_data[1][self.xy_coords[index_x][1], self.xy_coords[index_x][0]]

        return data

    #################
    # HYDRO
    # 3D
    #################
    def read_variable_sea_water_temperature_at_time_and_depth(self,index_t,index_z):
        data = np.zeros([self.nbPoints])
        data[:] = np.nan
        all_data = self.coverage.reader.read_variable_sea_water_temperature_at_time_and_depth(index_t, index_z,0,self.coverage.get_x_size(),0,self.coverage.get_y_size())

        for index_x in range(0, self.nbPoints):
            data[index_x] = all_data[self.xy_coords[index_x][1], self.xy_coords[index_x][0]]

        return data

    def read_variable_sea_water_salinity_at_time_and_depth(self,index_t,index_z):
        data = np.zeros([self.nbPoints])
        data[:] = np.nan

        all_data = self.coverage.reader.read_variable_sea_water_salinity_at_time_and_depth(index_t, index_z,0,self.coverage.get_x_size(),0,self.coverage.get_y_size())

        for index_x in range(0, self.nbPoints):
            data[index_x] = all_data[self.xy_coords[index_x][1], self.xy_coords[index_x][0]]

        return data

    def read_variable_baroclinic_sea_water_velocity_at_time_and_depth(self,index_t,index_z):
        data = np.zeros([2, self.nbPoints])
        data[:] = np.nan

        all_data = self.coverage.reader.read_variable_baroclinic_sea_water_velocity_at_time_and_depth(index_t, index_z)

        for index_x in range(0, self.nbPoints):
            data[0][index_x] = all_data[0][self.xy_coords[index_x][1], self.xy_coords[index_x][0]]
            data[1][index_x] = all_data[1][self.xy_coords[index_x][1], self.xy_coords[index_x][0]]

        return data

    #################
    # WAVES
    # Sea Surface
    #################
    def read_variable_sea_surface_wave_significant_height_at_time(self,index_t):
        data = np.zeros([self.nbPoints])
        data[:] = np.nan

        all_data = self.coverage.reader.read_variable_sea_surface_wave_significant_height_at_time(index_t)

        for index_x in range(0, self.nbPoints):
            data[index_x] = all_data[self.xy_coords[index_x][1], self.xy_coords[index_x][0]]

        return data

    def read_variable_sea_surface_wave_mean_period_at_time(self,index_t):
        data = np.zeros([self.nbPoints])
        data[:] = np.nan

        all_data = self.coverage.reader.read_variable_sea_surface_wave_mean_period_at_time(index_t)

        for index_x in range(0, self.nbPoints):
            data[index_x] = all_data[self.xy_coords[index_x][1], self.xy_coords[index_x][0]]

        return data

    def read_variable_sea_surface_wave_to_direction_at_time(self,index_t):
        data = np.zeros([self.nbPoints])
        data[:] = np.nan

        all_data = self.coverage.reader.read_variable_sea_surface_wave_to_direction_at_time(index_t)

        for index_x in range(0, self.nbPoints):
            data[index_x] = all_data[self.xy_coords[index_x][1], self.xy_coords[index_x][0]]

        return data

    def read_variable_sea_surface_wave_stokes_drift_velocity_at_time(self,index_t):
        data = np.zeros([2, self.nbPoints])
        data[:] = np.nan

        all_data = self.coverage.reader.read_variable_sea_surface_wave_stokes_drift_velocity_at_time(index_t)

        for index_x in range(0, self.nbPoints):
            data[0][index_x] = all_data[0][self.xy_coords[index_x][1], self.xy_coords[index_x][0]]
            data[1][index_x] = all_data[1][self.xy_coords[index_x][1], self.xy_coords[index_x][0]]

        return data

    #################
    # WAVES
    # Momentum flux
    #################
    def read_variable_atmosphere_momentum_flux_to_waves_at_time(self,index_t):
        data = np.zeros([2, self.nbPoints])
        data[:] = np.nan

        all_data = self.coverage.reader.read_variable_atmosphere_momentum_flux_to_waves_at_time(index_t)

        for index_x in range(0, self.nbPoints):
            data[0][index_x] = all_data[0][self.xy_coords[index_x][1], self.xy_coords[index_x][0]]
            data[1][index_x] = all_data[1][self.xy_coords[index_x][1], self.xy_coords[index_x][0]]

        return data

    def read_variable_waves_momentum_flux_to_ocean_at_time(self,index_t):
        data = np.zeros([2, self.nbPoints])
        data[:] = np.nan

        all_data = self.coverage.reader.read_variable_waves_momentum_flux_to_ocean_at_time(index_t)

        for index_x in range(0, self.nbPoints):
            data[0][index_x] = all_data[0][self.xy_coords[index_x][1], self.xy_coords[index_x][0]]
            data[1][index_x] = all_data[1][self.xy_coords[index_x][1], self.xy_coords[index_x][0]]

        return data

    #################
    # METEO
    # Surface air
    #################
    def read_variable_wind_stress_at_time(self,index_t):
        data = np.zeros([2, self.nbPoints])
        data[:] = np.nan

        all_data = self.coverage.reader.read_variable_wind_stress_at_time(index_t)

        for index_x in range(0, self.nbPoints):
            data[0][index_x] = all_data[0][self.xy_coords[index_x][1], self.xy_coords[index_x][0]]
            data[1][index_x] = all_data[1][self.xy_coords[index_x][1], self.xy_coords[index_x][0]]

        return data

    #################
    # METEO
    # At 10 m
    #################
    def read_variable_wind_10m_at_time(self,index_t):
        data = np.zeros([2, self.nbPoints])
        data[:] = np.nan

        all_data = self.coverage.reader.read_variable_wind_10m_at_time(index_t)

        for index_x in range(0, self.nbPoints):
            data[0][index_x] = all_data[0][self.xy_coords[index_x][1], self.xy_coords[index_x][0]]
            data[1][index_x] = all_data[1][self.xy_coords[index_x][1], self.xy_coords[index_x][0]]

        return data






