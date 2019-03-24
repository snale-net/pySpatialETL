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
from numpy import float32,float64,int32

import numpy as np
import logging

class DefaultMultiPointWriter(MultiPointWriter):

    def __init__(self, p,myFile):
        MultiPointWriter.__init__(self,p,myFile)

        if not isinstance(self.points, MultiPoint):
            raise ValueError("This writer supports only a MultiPoint object")

        self.ncfile = Dataset(self.filename, 'w', format='NETCDF4')
        self.ncfile.description = 'Generated with pyGeoSpatialETL'

        self.ncfile.data_source = str(self.points.data_source)
        self.ncfile.meta_data = str(self.points.meta_data)

        # Geo-points dimension
        self.ncfile.createDimension('point', self.points.get_nb_points())
        var = self.ncfile.createVariable('point', int32, ('point',))
        var.long_name = "point number";
        var.standard_name = "point_number";
        var.axis = "X";
        var[:] = range(0,self.points.get_nb_points());

        var = self.ncfile.createVariable('longitude', float32, ('point',), fill_value=9.96921e+36)
        var.long_name = "longitude";
        var.standard_name = "longitude";
        var.valid_min = "-180.0";
        var.valid_max = "180.0";
        var.units = VariableUnits.CANONICAL_UNITS[var.standard_name];
        var[:] = self.points.read_axis_x()

        var = self.ncfile.createVariable('latitude', float32, ('point',), fill_value=9.96921e+36)
        var.long_name = "latitude";
        var.standard_name = "latitude";
        var.valid_min = "-90.0";
        var.valid_max = "90.0";
        var.units = VariableUnits.CANONICAL_UNITS[var.standard_name];
        var[:] = self.points.read_axis_y()

    def close(self):
        self.ncfile.close()

    def write_variable_time(self):
        times = self.ncfile.createVariable('point_time', float64, ('point',))
        times.units = 'seconds since 1970-01-01 00:00:00'
        times.calendar = 'gregorian'
        times.standard_name = 'time'
        times.conventions = "UTC time"

        times[:] = date2num(self.points.read_variable_time(), units=times.units, calendar=times.calendar)

    def write_variable_bathymetry(self):
        var = self.ncfile.createVariable('bathymetry', float32, ('point',),
                                         fill_value=9.96921e+36)
        var.long_name = "Bathymetry";
        var.standard_name = "bathymetry";
        var.units = VariableUnits.CANONICAL_UNITS[var.standard_name];
        var.positive = "down" ;

        logging.info('[DefaultMultiPointWriter] Writing variable \'' + var.long_name + '\'')
        var[:] = self.points.read_variable_bathymetry()

    def write_variable_sea_surface_temperature(self):
        var = self.ncfile.createVariable('sea_surface_temperature', float32, ('point',),
                                         fill_value=9.96921e+36)
        var.long_name = "Sea Surface Temperature";
        var.standard_name = "sea_surface_temperature";
        var.units = VariableUnits.CANONICAL_UNITS[var.standard_name];

        logging.info('[DefaultMultiPointWriter] Writing variable \'' + var.long_name + '\'')
        var[:] = self.points.read_variable_sea_surface_temperature()

    def write_variable_sea_water_electrical_conductivity(self):
        var = self.ncfile.createVariable('sea_water_electrical_conductivity', float32, ('point',),
                                         fill_value=9.96921e+36)
        var.long_name = "Sea Water Electrical Conductivity";
        var.standard_name = "sea_water_electrical_conductivity";
        var.units = VariableUnits.CANONICAL_UNITS[var.standard_name];

        logging.info('[DefaultMultiPointWriter] Writing variable \'' + var.long_name + '\'')
        var[:] = self.points.read_variable_sea_water_electrical_conductivity()

    def write_variable_sea_water_pressure_at_sea_water_surface(self):
        var = self.ncfile.createVariable('sea_water_pressure_at_sea_water_surface', float32, ('point'),
                                         fill_value=9.96921e+36)
        var.long_name = "Sea Water Pressure At Sea Water Surface";
        var.standard_name = "sea_water_pressure_at_sea_water_surface";
        var.units = VariableUnits.CANONICAL_UNITS[var.standard_name];

        logging.info('[DefaultMultiPointWriter] Writing variable \'' + var.long_name + '\'')
        var[:] = self.points.read_variable_sea_water_pressure_at_sea_water_surface()

    def write_variable_sea_surface_salinity(self):
        var = self.ncfile.createVariable('sea_surface_salinity', float32, ('point'), fill_value=9.96921e+36)
        var.long_name = "Sea Surface Salinity";
        var.standard_name = "sea_surface_salinity";
        var.units = VariableUnits.CANONICAL_UNITS[var.standard_name];

        logging.info('[DefaultMultiPointWriter] Writing variable \'' + var.long_name + '\'')
        var[:] = self.points.read_variable_sea_surface_salinity()

    def write_variable_sea_surface_density(self):
        var = self.ncfile.createVariable('sea_surface_density', float32, ('point',), fill_value=9.96921e+36)
        var.long_name = "Sea Surface Density";
        var.standard_name = "sea_surface_densitiy";
        var.units = VariableUnits.CANONICAL_UNITS[var.standard_name];

        logging.info('[DefaultMultiPointWriter] Writing variable \'' + var.long_name + '\'')
        var[:] = self.points.read_variable_sea_surface_density()

    def write_variable_sea_water_turbidity(self):
        var = self.ncfile.createVariable('sea_water_turbidity', float32, ('point',), fill_value=9.96921e+36)
        var.long_name = "Sea Surface Turbidity";
        var.standard_name = "sea_water_turbidity";
        var.units = VariableUnits.CANONICAL_UNITS[var.standard_name];

        logging.info('[DefaultMultiPointWriter] Writing variable \'' + var.long_name + '\'')
        var[:] = self.points.read_variable_sea_water_turbidity()
