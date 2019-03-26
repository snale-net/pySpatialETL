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
from utils.VariableDefinition import VariableDefinition
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
        self.ncfile.createDimension(VariableDefinition.VARIABLE_NAME['point'], self.points.get_nb_points())
        var = self.ncfile.createVariable(VariableDefinition.VARIABLE_NAME['point'], int32,
                                         (VariableDefinition.VARIABLE_NAME['point'],))
        var.long_name = VariableDefinition.LONG_NAME['point']
        var.standard_name = VariableDefinition.STANDARD_NAME['point']
        var.axis = "X";
        var.units = VariableDefinition.CANONICAL_UNITS['point'];
        var[:] = range(0, self.points.get_nb_points());

        var = self.ncfile.createVariable(VariableDefinition.VARIABLE_NAME['latitude'], float32,
                                         (VariableDefinition.VARIABLE_NAME['point'],), fill_value=9.96921e+36)
        var.long_name = VariableDefinition.LONG_NAME['latitude']
        var.standard_name = VariableDefinition.STANDARD_NAME['latitude']
        var.valid_min = "-90.0";
        var.valid_max = "90.0";
        var.units = VariableDefinition.CANONICAL_UNITS['latitude']
        var[:] = self.points.read_axis_y()

        var = self.ncfile.createVariable(VariableDefinition.VARIABLE_NAME['longitude'], float32,
                                         (VariableDefinition.VARIABLE_NAME['point'],), fill_value=9.96921e+36)
        var.long_name = VariableDefinition.LONG_NAME['longitude']
        var.standard_name = VariableDefinition.STANDARD_NAME['longitude']
        var.valid_min = "-180.0";
        var.valid_max = "180.0";
        var.units = VariableDefinition.CANONICAL_UNITS['longitude']
        var[:] = self.points.read_axis_x()

    def close(self):
        self.ncfile.close()

    def write_variable_time(self):
        times = self.ncfile.createVariable(VariableDefinition.VARIABLE_NAME['time'], float64, (VariableDefinition.VARIABLE_NAME['point'],))
        times.units = 'seconds since 1970-01-01 00:00:00'
        times.calendar = 'gregorian'
        times.standard_name = 'time'
        times.conventions = "UTC time"

        times[:] = date2num(self.points.read_variable_time(), units=times.units, calendar=times.calendar)

    def write_variable_bathymetry(self):
        var = self.ncfile.createVariable(VariableDefinition.VARIABLE_NAME['bathymetry'], float32, (VariableDefinition.VARIABLE_NAME['point'],),
                                         fill_value=9.96921e+36)
        var.long_name = VariableDefinition.LONG_NAME['bathymetry']
        var.standard_name = VariableDefinition.STANDARD_NAME['bathymetry']
        var.units = VariableDefinition.CANONICAL_UNITS['bathymetry']
        var.positive = "down" ;

        logging.info('[DefaultMultiPointWriter] Writing variable \'' + str(VariableDefinition.LONG_NAME['bathymetry']) + '\'')
        var[:] = self.points.read_variable_bathymetry()

    def write_variable_sea_surface_temperature(self):
        var = self.ncfile.createVariable(VariableDefinition.VARIABLE_NAME['sea_surface_temperature'], float32, (VariableDefinition.VARIABLE_NAME['point'],),
                                         fill_value=9.96921e+36)
        var.long_name = VariableDefinition.LONG_NAME['sea_surface_temperature']
        var.standard_name = VariableDefinition.STANDARD_NAME['sea_surface_temperature']
        var.units = VariableDefinition.CANONICAL_UNITS['sea_surface_temperature']

        logging.info('[DefaultMultiPointWriter] Writing variable \'' + str(VariableDefinition.LONG_NAME['sea_surface_temperature']) + '\'')
        var[:] = self.points.read_variable_sea_surface_temperature()

    def write_variable_sea_water_electrical_conductivity(self):
        var = self.ncfile.createVariable(VariableDefinition.VARIABLE_NAME['sea_water_electrical_conductivity'], float32, (VariableDefinition.VARIABLE_NAME['point'],),
                                         fill_value=9.96921e+36)
        var.long_name = VariableDefinition.LONG_NAME['sea_water_electrical_conductivity']
        var.standard_name = VariableDefinition.STANDARD_NAME['sea_water_electrical_conductivity']
        var.units = VariableDefinition.CANONICAL_UNITS['sea_water_electrical_conductivity']

        logging.info('[DefaultMultiPointWriter] Writing variable \'' + str(VariableDefinition.LONG_NAME['sea_water_electrical_conductivity']) + '\'')
        var[:] = self.points.read_variable_sea_water_electrical_conductivity()

    def write_variable_sea_water_pressure_at_sea_water_surface(self):
        var = self.ncfile.createVariable(VariableDefinition.VARIABLE_NAME['sea_water_pressure_at_sea_water_surface'], float32, (VariableDefinition.VARIABLE_NAME['point']),
                                         fill_value=9.96921e+36)
        var.long_name = VariableDefinition.LONG_NAME['sea_water_pressure_at_sea_water_surface']
        var.standard_name = VariableDefinition.STANDARD_NAME['sea_water_pressure_at_sea_water_surface']
        var.units = VariableDefinition.CANONICAL_UNITS['sea_water_pressure_at_sea_water_surface']

        logging.info('[DefaultMultiPointWriter] Writing variable \'' + str( VariableDefinition.LONG_NAME['sea_water_pressure_at_sea_water_surface']) + '\'')
        var[:] = self.points.read_variable_sea_water_pressure_at_sea_water_surface()

    def write_variable_sea_surface_salinity(self):
        var = self.ncfile.createVariable(VariableDefinition.VARIABLE_NAME['sea_surface_salinity'], float32, (VariableDefinition.VARIABLE_NAME['point']), fill_value=9.96921e+36)
        var.long_name = VariableDefinition.LONG_NAME['sea_surface_salinity']
        var.standard_name = VariableDefinition.STANDARD_NAME['sea_surface_salinity']
        var.units = VariableDefinition.CANONICAL_UNITS['sea_surface_salinity']

        logging.info('[DefaultMultiPointWriter] Writing variable \'' + str(VariableDefinition.LONG_NAME['sea_surface_salinity']) + '\'')
        var[:] = self.points.read_variable_sea_surface_salinity()

    def write_variable_sea_surface_density(self):
        var = self.ncfile.createVariable(VariableDefinition.VARIABLE_NAME['sea_surface_density'], float32, (VariableDefinition.VARIABLE_NAME['point'],), fill_value=9.96921e+36)
        var.long_name = VariableDefinition.LONG_NAME['sea_surface_density']
        var.standard_name = VariableDefinition.STANDARD_NAME['sea_surface_density']
        var.units =VariableDefinition.CANONICAL_UNITS['sea_surface_density']

        logging.info('[DefaultMultiPointWriter] Writing variable \'' + str(VariableDefinition.LONG_NAME['sea_surface_density']) + '\'')
        var[:] = self.points.read_variable_sea_surface_density()

    def write_variable_sea_water_turbidity(self):
        var = self.ncfile.createVariable(VariableDefinition.VARIABLE_NAME['sea_water_turbidity'], float32, (VariableDefinition.VARIABLE_NAME['point'],), fill_value=9.96921e+36)
        var.long_name = VariableDefinition.LONG_NAME['sea_water_turbidity']
        var.standard_name = VariableDefinition.STANDARD_NAME['sea_water_turbidity']
        var.units = VariableDefinition.CANONICAL_UNITS['sea_water_turbidity']

        logging.info('[DefaultMultiPointWriter] Writing variable \'' + str(VariableDefinition.LONG_NAME['sea_water_turbidity']) + '\'')
        var[:] = self.points.read_variable_sea_water_turbidity()
