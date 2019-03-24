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
from point.io.MultiPointWriter import MultiPointWriter
from point.MultiPoint import MultiPoint
from utils.VariableUnits import VariableUnits
from netCDF4 import Dataset
from netCDF4 import date2num
from numpy import float32
from numpy import float64
import numpy as np
import logging
import pandas


class DefaultMultiPointWriter(MultiPointWriter):

    def __init__(self, s,myFile):
        MultiPointWriter.__init__(self, s, myFile);

        if not isinstance(self.points, MultiPoint):
            raise ValueError("This writer supports only a MultiPoint object")

        index = pandas.Index(range(0,self.points.get_nb_points()))
        self.data = pandas.DataFrame(index=index)

        self.data['longitude'] = self.points.read_axis_x()
        self.data['latitude'] = self.points.read_axis_y()

    def close(self):
        self.data.to_csv(self.filename, sep='\t', columns=list(self.data), header=False, encoding='utf-8', na_rep="NaN")

        file = open(self.filename, "r+")
        old = file.read()  # read everything in the file
        file.seek(0)  # rewind

        file.write("############################################################ \n\
# Coordinate Reference System : WGS84 \n\
# Data source : " + str(self.points.data_source) + " \n\
# Meta Data : " + str(self.points.meta_data) + " \n\
# Separator: Tabulation \\t \n\
# Column 1: Point number \n")

        column = 2
        for key in list(self.data):
            file.write("# Column " + str(column) + ": " + str(key) + " " + str(
                VariableUnits.CANONICAL_UNITS[key]) + " FillValue: NaN \n")
            column = column + 1

        file.write("# Generated with pyGeoSpatialETL\n")
        file.write("############################################################\n")

        file.write(old)  # write the new line before
        file.close()

    def write_variable_sea_surface_height(self):
        logging.info('[DefaultMultiPointWriter] Writing variable \'Sea Surface Height\'')
        self.data['sea_surface_height'] = self.points.read_variable_sea_surface_height()

    def write_variable_sea_surface_temperature(self):
        logging.info('[DefaultMultiPointWriter] Writing variable \'Sea Surface Temperature\'')
        self.data['sea_surface_temperature'] = self.points.read_variable_sea_surface_temperature()

    def write_variable_sea_surface_wave_significant_height(self):
        logging.info('[DefaultMultiPointWriter] Writing variable \'Sea Surface Wave Significant Height\'')
        self.data['sea_surface_wave_significant_height'] = self.points.read_variable_sea_surface_wave_significant_height()

    def write_variable_sea_surface_wave_mean_period(self):
        logging.info('[DefaultMultiPointWriter] Writing variable \'Sea Surface Wave Mean Period\'')
        self.data['sea_surface_wave_mean_periode'] = self.points.read_variable_sea_surface_wave_mean_period()

    def write_variable_sea_surface_salinity(self):
        logging.info('[DefaultMultiPointWriter] Writing variable \'Sea Surface Salinity\'')
        self.data['sea_surface_salinity'] = self.points.read_variable_sea_surface_salinity()

    def write_variable_sea_surface_density(self):
        logging.info('[DefaultMultiPointWriter] Writing variable \'Sea Surface Density\'')
        self.data['sea_surface_density'] = self.points.read_variable_sea_surface_density()

    def write_variable_sea_water_turbidity(self):
        logging.info('[DefaultMultiPointWriter] Writing variable \'Sea Water Turbidity\'')
        self.data['sea_water_turbidity'] = self.points.read_variable_sea_water_turbidity()

    def write_variable_bathymetry(self):
        logging.info('[DefaultMultiPointWriter] Writing variable \'Bathymetry\'')
        self.data['bathymetry'] = self.points.read_variable_bathymetry()

    def write_variable_wind_speed_10m(self):
        logging.info('[DefaultMultiPointWriter] Writing variable \'Wind Speed 10m\'')
        self.data['wind_speed_10m'] = self.points.read_variable_wind_speed_10m()

    def write_variable_wind_from_direction_10m(self):
        logging.info('[DefaultMultiPointWriter] Writing variable \'Wind From Direction 10m\'')
        self.data['wind_from_direction_10m'] = self.points.read_variable_wind_from_direction_10m()

    def write_variable_surface_air_pressure(self):
        logging.info('[DefaultMultiPointWriter] Writing variable \'Surface Air Pressure\'')
        self.data['surface_air_pressure'] = self.points.read_variable_surface_air_pressure()

    def write_variable_rainfall_amount(self):
        logging.info('[DefaultMultiPointWriter] Writing variable \'Rainfall Amount\'')
        self.data['rainfall_amount'] = self.points.read_variable_rainfall_amount()

    def write_variable_sea_water_pressure_at_sea_water_surface(self):
        logging.info('[DefaultMultiPointWriter] Writing variable \Sea Water Pressure At Sea Water Surface\'')
        self.data['sea_water_pressure_at_sea_water_surface'] = self.points.read_variable_sea_water_pressure_at_sea_water_surface()

    def write_variable_sea_water_electrical_conductivity(self):
        logging.info('[DefaultMultiPointWriter] Writing variable \'Sea Water Electrical Conductivity\'')
        self.data['sea_water_electrical_conductivity'] = self.points.read_variable_sea_water_electrical_conductivity()
