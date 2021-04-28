##! /usr/bin/env python2.7
# -*- coding: utf-8 -*-
#
# pySpatialETL is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# any later version.
#
# pySpatialETL is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# Author : Fabien Rétif - fabien.retif@zoho.com
#
from __future__ import division, print_function, absolute_import

import glob
import os
import re

import cftime
import numpy as np
from osgeo import gdal

from spatialetl.coverage.TimeCoverage import TimeCoverage
from spatialetl.coverage.io.CoverageReader import CoverageReader
from spatialetl.utils.path import path_leaf


class INSPIREReader(CoverageReader):

    VARIABLES = [
        'BRIGHTNESS_TEMPERATURE__GROUND_OR_WATER_SURFACE',
        'TOTAL_WATER_PRECIPITATION__GROUND_OR_WATER_SURFACE',
        'DOWNWARD_SHORT_WAVE_RADIATION_FLUX__GROUND_OR_WATER_SURFACE',
        'SHORT_WAVE_RADIATION_FLUX__GROUND_OR_WATER_SURFACE',
        'RELATIVE_HUMIDITY__SPECIFIC_HEIGHT_LEVEL_ABOVE_GROUND',
        'TOTAL_SNOW_PRECIPITATION__GROUND_OR_WATER_SURFACE',
        'PRESSURE__MEAN_SEA_LEVEL',
        'DEW_POINT_TEMPERATURE__SPECIFIC_HEIGHT_LEVEL_ABOVE_GROUND',
        'TEMPERATURE__GROUND_OR_WATER_SURFACE',
        'U_COMPONENT_OF_WIND__SPECIFIC_HEIGHT_LEVEL_ABOVE_GROUND',
        'V_COMPONENT_OF_WIND__SPECIFIC_HEIGHT_LEVEL_ABOVE_GROUND']

    VARIABLES_PT3H = [
        'TOTAL_WATER_PRECIPITATION__GROUND_OR_WATER_SURFACE',
        'DOWNWARD_SHORT_WAVE_RADIATION_FLUX__GROUND_OR_WATER_SURFACE',
        'SHORT_WAVE_RADIATION_FLUX__GROUND_OR_WATER_SURFACE',
        'TOTAL_SNOW_PRECIPITATION__GROUND_OR_WATER_SURFACE']

    def __init__(self, myFile):
        CoverageReader.__init__(self,myFile);
        if os.path.isfile(self.filename):
            self.files= [self.filename]
        elif os.path.isdir(self.filename):
            self.files = sorted(glob.glob(os.path.join(self.filename, "*.tiff")))
        elif self.filename.endswith("*"):
            self.files = sorted(glob.glob(self.filename + ".tiff"))
        else:
            raise ValueError("Unable to decode file " + str(self.filename))

        self.times = []
        self.variable_files = {}
        for var in INSPIREReader.VARIABLES:
            self.variable_files[var] = []

        for file in self.files:
            groups = re.search("^([0-9]{4})([0-9]{2})([0-9]{2})\_([0-9]{2})([0-9]{2})([0-9]{2})\_MF\_([A-Z_]+).tiff$",path_leaf(file))
            if groups:
                current_time = cftime.datetime(int(groups.group(1)), int(groups.group(2)), int(groups.group(3)), int(groups.group(4)), int(groups.group(5)),
                                    int(groups.group(6)))
                if current_time not in self.times:
                    self.times.append(current_time)

                if self.variable_files[groups.group(7)] is None:
                    self.variable_files[groups.group(7)] = []

                self.variable_files[groups.group(7)].append(file)

        if len(self.times) == 0:
            raise ValueError("Unable to find Inspire raw output filename.")

        if len(self.times)*len(INSPIREReader.VARIABLES) != len(self.files):
            print("Nb times classiques : "+str(len(self.times)))
            print("Nb variables : " + str(len(INSPIREReader.VARIABLES)))
            print("Nb total times : " + str(len(self.times) * len(INSPIREReader.VARIABLES)))
            print("Nb files : " + str(len(self.files)))
            raise ValueError("Tous les variables n'ont pas les mêmes temps :")

        self.last_opened_t_index = 0

        self.tifffile = gdal.Open(self.variable_files["PRESSURE__MEAN_SEA_LEVEL"][self.last_opened_t_index])
        self.y = None
        self.x = None

    def open_file(self,varname, index_t):
        self.close()
        self.tifffile = gdal.Open(self.variable_files[varname][index_t])

    def close(self):
        self.tifffile = None

    def is_regular_grid(self):
        return True

    def pixel2coord(self, y, x):
        # unravel GDAL affine transform parameters
        c, a, b, f, d, e = self.tifffile.GetGeoTransform()
        """Returns global coordinates to pixel center using base-0 raster index"""
        xp = a * y + b * x + a * 0.5 + b * 0.5 + c
        yp = d * y + e * x + d * 0.5 + e * 0.5 + f
        return (xp, yp)

    def get_x_size(self):
        return self.tifffile.RasterXSize

    def get_y_size(self):
        return self.tifffile.RasterYSize

    def get_t_size(self):
        return len(self.times)

    # Axis
    def read_axis_t(self,tmin,tmax,timestamp):

        if timestamp == 1:
            return [(t - TimeCoverage.TIME_DATUM).total_seconds() \
                    for t in self.times[tmin:tmax]];
        else:
            return self.times[tmin:tmax]
        
    def read_axis_x(self,xmin,xmax,ymin,ymax):

        if self.x is None:
            width = self.tifffile.RasterXSize

            self.x = np.zeros([width])
            for x in range(0,width):
                self.x[x] = self.pixel2coord(x,0)[0]

        return self.x[xmin:xmax]

    def read_axis_y(self,xmin,xmax,ymin,ymax):
        if self.y is None:
            height = self.tifffile.RasterYSize

            self.y = np.zeros([height])
            for y in range(0,height):
                self.y[y] = self.pixel2coord(0,y)[1]

        return self.y[ymin:ymax]

    # Scalar
    def read_variable_2D_sea_binary_mask(self, xmin, xmax, ymin, ymax):
        data = np.zeros([ymax-ymin,xmax-xmin])
        data[:] = 1
        return data

    def read_variable_2D_sea_binary_mask_at_time(self,intex_t, xmin, xmax, ymin, ymax):
        data = np.zeros([ymax - ymin, xmax - xmin])
        data[:] = 1
        return data

    def read_variable_2D_land_binary_mask(self, xmin, xmax, ymin, ymax):
        data = np.zeros([ymax-ymin,xmax-xmin])
        data[:] = 1
        return data

    def read_variable_2D_land_binary_mask_at_time(self,intex_t, xmin, xmax, ymin, ymax):
        data = np.zeros([ymax - ymin, xmax - xmin])
        data[:] = 1
        return data

    #################
    # METEO
    # 2D
    #################

    def read_variable_rainfall_amount_at_time(self, index_t, xmin, xmax, ymin, ymax):
        self.open_file("TOTAL_WATER_PRECIPITATION__GROUND_OR_WATER_SURFACE",index_t)
        band = self.tifffile.GetRasterBand(1)
        data = band.ReadAsArray()
        return np.ma.filled(data[ymin:ymax, xmin:xmax], fill_value=np.nan)

    #################
    # METEO
    # Sea surface
    #################

    def read_variable_surface_air_pressure_at_time(self, index_t, xmin, xmax, ymin, ymax):
        self.open_file("PRESSURE__MEAN_SEA_LEVEL", index_t)
        band = self.tifffile.GetRasterBand(1)
        data = band.ReadAsArray()
        data *= 0.01  # Pa to hPa
        return np.ma.filled(data[ymin:ymax, xmin:xmax], fill_value=np.nan)



    def read_variable_sea_surface_air_pressure_at_time(self, index_t, xmin, xmax, ymin, ymax):
        self.open_file("PRESSURE__MEAN_SEA_LEVEL", index_t)
        band = self.tifffile.GetRasterBand(1)
        data = band.ReadAsArray()
        data *= 0.01  # Pa to hPa
        return np.ma.filled(data[ymin:ymax, xmin:xmax], fill_value=np.nan)

    def read_variable_surface_downward_sensible_heat_flux_at_time(self, index_t, xmin, xmax, ymin, ymax):
        data = np.zeros([ymax - ymin, xmax - xmin])
        data[:] = -50530
        return data

    def read_variable_surface_downward_latent_heat_flux_at_time(self, index_t, xmin, xmax, ymin, ymax):
        data = np.zeros([ymax - ymin, xmax - xmin])
        data[:] = -880170
        return data

    def read_variable_surface_air_temperature_at_time(self, index_t, xmin, xmax, ymin, ymax):
        self.open_file("TEMPERATURE__GROUND_OR_WATER_SURFACE", index_t)
        band = self.tifffile.GetRasterBand(1)
        data = band.ReadAsArray()
        #data -= 273.15 # Kelvin to Celsius
        #data = convert_temperature(band.ReadAsArray(), 'Kelvin', 'Celsius')
        return np.ma.filled(data[ymin:ymax, xmin:xmax], fill_value=np.nan)

    def read_variable_dew_point_temperature_at_time(self, index_t, xmin, xmax, ymin, ymax):
        self.open_file("DEW_POINT_TEMPERATURE__SPECIFIC_HEIGHT_LEVEL_ABOVE_GROUND", index_t)
        band = self.tifffile.GetRasterBand(1)
        data = band.ReadAsArray()
        #data -= 273.15 # Kelvin to Celsius
        #data = convert_temperature(band.ReadAsArray(), 'Kelvin','Celsius')
        return np.ma.filled(data[ymin:ymax, xmin:xmax], fill_value=np.nan)

    def read_variable_surface_downward_solar_radiation_at_time(self, index_t, xmin, xmax, ymin, ymax):
        data = np.zeros([ymax - ymin, xmax - xmin])
        data[:] = 4.19e+6
        return data

    def read_variable_surface_downward_thermal_radiation_at_time(self, index_t, xmin, xmax, ymin, ymax):
        data = np.zeros([ymax - ymin, xmax - xmin])
        data[:] = 4.67e+6
        return data

    def read_variable_surface_solar_radiation_at_time(self, index_t, xmin, xmax, ymin, ymax):
        data = np.zeros([ymax - ymin, xmax - xmin])
        data[:] = 4e+6
        return data

    def read_variable_surface_thermal_radiation_at_time(self, index_t, xmin, xmax, ymin, ymax):
        data = np.zeros([ymax - ymin, xmax - xmin])
        data[:] = -303395
        return data

    #################
    # METEO
    # At 10 m
    #################

    def read_variable_wind_10m_at_time(self,index_t,xmin,xmax,ymin,ymax):

        self.open_file("U_COMPONENT_OF_WIND__SPECIFIC_HEIGHT_LEVEL_ABOVE_GROUND", index_t)
        band = self.tifffile.GetRasterBand(1)
        u_comp = band.ReadAsArray()

        self.open_file("V_COMPONENT_OF_WIND__SPECIFIC_HEIGHT_LEVEL_ABOVE_GROUND", index_t)
        band = self.tifffile.GetRasterBand(1)
        v_comp = band.ReadAsArray()

        return [u_comp[ymin:ymax,xmin:xmax],v_comp[ymin:ymax,xmin:xmax]]


