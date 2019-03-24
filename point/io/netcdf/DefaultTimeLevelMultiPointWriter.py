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
from point.TimeLevelMultiPoint import TimeLevelMultiPoint
from utils.VariableUnits import VariableUnits
from netCDF4 import Dataset
from netCDF4 import date2num
from numpy import float32,float64,int32

import numpy as np
import logging

class DefaultTimeLevelMultiPointWriter(MultiPointWriter):

    def __init__(self, p,myFile):
        MultiPointWriter.__init__(self,p,myFile)

        if not isinstance(self.points, TimeLevelMultiPoint):
            raise ValueError("This writer supports only a TimeLevelMultiPoint object")

        self.ncfile = Dataset(self.filename, 'w', format='NETCDF4')
        self.ncfile.description = 'Generated with pyGeoSpatialETL'

        self.ncfile.data_source = str(self.points.data_source)
        self.ncfile.meta_data = str(self.points.meta_data)

        # Profil dimension
        self.ncfile.createDimension('point', self.points.get_nb_points())
        var = self.ncfile.createVariable('point', int32, ('point',))
        var.long_name = "Point number";
        var.standard_name = "point_number";
        var.axis = "X";
        var[:] = range(0, self.points.get_nb_points());

        var = self.ncfile.createVariable('latitude', float32, ('point',), fill_value=9.96921e+36)
        var.long_name = "latitude";
        var.standard_name = "latitude";
        var.valid_min = "-90.0";
        var.valid_max = "90.0";
        var.units = VariableUnits.CANONICAL_UNITS[var.standard_name];
        var[:] = self.points.read_axis_y()

        var = self.ncfile.createVariable('longitude', float32, ('point',), fill_value=9.96921e+36)
        var.long_name = "longitude";
        var.standard_name = "longitude";
        var.valid_min = "-180.0";
        var.valid_max = "180.0";
        var.units = VariableUnits.CANONICAL_UNITS[var.standard_name];
        var[:] = self.points.read_axis_x()

        # Time dimension
        self.ncfile.createDimension('time', None)
        times = self.ncfile.createVariable('time', float64, ('time',))
        times.units= 'seconds since 1970-01-01 00:00:00'
        times.calendar= 'gregorian'
        times.standard_name= 'time'
        times.axis='T'
        times.conventions = "UTC time"

        times[:] = date2num(self.points.read_axis_t(), units = times.units, calendar = times.calendar)

        # Depth dimension
        self.ncfile.createDimension('depth', self.points.get_z_size())
        var = self.ncfile.createVariable('depth', float64, ('depth',))
        var.standard_name = 'depth'
        var.long_name = "Depth"
        var.positive = "down";
        var.axis = 'Z'
        var.units = VariableUnits.CANONICAL_UNITS[var.standard_name];
        var[:] = self.points.read_axis_z()


    def close(self):
        self.ncfile.close()

    def write_variable_baroclinic_sea_water_velocity(self):

        ucur = self.ncfile.createVariable('baroclinic_eastward_sea_water_velocity', float32, ('time','depth','point',),
                                          fill_value=9.96921e+36)
        ucur.long_name = "Baroclinic Eastward Sea Water Velocity";
        ucur.standard_name = "baroclinic_eastward_sea_water_velocity";
        ucur.units = VariableUnits.CANONICAL_UNITS[ucur.standard_name];
        ucur.comment = "cur=sqrt(U**2+V**2)";

        vcur = self.ncfile.createVariable('baroclinic_northward_sea_water_velocity', float32, ('time','depth','point'),
                                          fill_value=9.96921e+36)
        vcur.long_name = "Baroclinic Northward Sea Water Velocity";
        vcur.standard_name = "baroclinic_northward_sea_water_velocity";
        vcur.units = VariableUnits.CANONICAL_UNITS[vcur.standard_name];
        vcur.comment = "cur=sqrt(U**2+V**2)";

        logging.info('[DefaultTimeLevelMultiPointWriter] Writing variable \'Baroclinic Sea Water Velocity\'\'')

        time_index = 0
        for time in self.points.read_axis_t():
            logging.debug(
                '[DefaultTimeLevelMultiPointWriter] Writing variable \'Baroclinic Sea Water Velocity\' at time \'' + str(time) + '\'')

            z_index = 0
            for depth in self.points.read_axis_z():

                cur = self.points.read_variable_baroclinic_sea_water_velocity_at_time_and_depth(time, depth)

                ucur[time_index:time_index + 1,z_index:z_index + 1] = cur[0]
                vcur[time_index:time_index + 1,z_index:z_index + 1] = cur[1]

                z_index = z_index + 1

            time_index += 1


    def write_variable_sea_water_temperature(self):

        var = self.ncfile.createVariable('sea_water_temperature', float32, ('time','depth', 'point'), fill_value=9.96921e+36)
        var.long_name = "Sea Water Temperature";
        var.standard_name = "sea_water_temperature";
        var.units = VariableUnits.CANONICAL_UNITS[var.standard_name];

        logging.info('[DefaultTimeLevelMultiPointWriter] Writing variable \'' + var.long_name + '\'')

        time_index = 0
        for time in self.points.read_axis_t():
            logging.debug('[DefaultTimeLevelMultiPointWriter] Writing variable \''+var.long_name+'\' at time \'' + str(time) + '\'')

            z_index = 0
            for depth in self.points.read_axis_z():
                data = self.points.read_variable_sea_water_temperature_at_time_and_depth(time,depth)

                var[time_index:time_index + 1,z_index:z_index + 1] = data
                z_index = z_index +1

            time_index += 1

    def write_variable_sea_water_salinity(self):
        var = self.ncfile.createVariable('sea_water_salinity', float32, ('time','depth', 'point'), fill_value=9.96921e+36)
        var.long_name = "Sea Water Salinity";
        var.standard_name = "sea_water_salinity";
        var.units = VariableUnits.CANONICAL_UNITS[var.standard_name];

        logging.info('[DefaultTimeLevelMultiPointWriter] Writing variable \'' + var.long_name + '\'')

        time_index = 0
        for time in self.points.read_axis_t():
            logging.debug(
                '[DefaultTimeLevelMultiPointWriter] Writing variable \'' + var.long_name + '\' at time \'' + str(time) + '\'')

            z_index = 0
            for depth in self.points.read_axis_z():
                data = self.points.read_variable_sea_water_salinity_at_time_and_depth(time, depth)

                var[time_index:time_index + 1, z_index:z_index + 1] = data
                z_index = z_index + 1

            time_index += 1

    def write_variable_sea_water_density(self):
        var = self.ncfile.createVariable('sea_water_density', float32, ('time','depth', 'point',), fill_value=9.96921e+36)
        var.long_name = "Sea Water Density";
        var.standard_name = "sea_water_densitiy";
        var.units = VariableUnits.CANONICAL_UNITS[var.standard_name];

        logging.info('[DefaultTimeLevelMultiPointWriter] Writing variable \'' + var.long_name + '\'')

        time_index = 0
        for time in self.points.read_axis_t():
            logging.debug(
                '[DefaultTimeLevelMultiPointWriter] Writing variable \'' + var.long_name + '\' at time \'' + str(time) + '\'')

            z_index = 0
            for depth in self.points.read_axis_z():
                data = self.points.read_variable_sea_water_density_at_time_and_depth(time, depth)

                var[time_index:time_index + 1, z_index:z_index + 1] = data
                z_index = z_index + 1

            time_index += 1

    def write_variable_sea_water_turbidity(self):
        var = self.ncfile.createVariable('sea_water_turbidity', float32, ('time','depth', 'point',), fill_value=9.96921e+36)
        var.long_name = "Sea Water Turbidity";
        var.standard_name = "sea_water_turbidity";
        var.units = VariableUnits.CANONICAL_UNITS[var.standard_name];

        logging.info('[DefaultTimeLevelMultiPointWriter] Writing variable \'' + var.long_name + '\'')

        time_index = 0
        for time in self.points.read_axis_t():
            logging.debug(
                '[DefaultTimeLevelMultiPointWriter] Writing variable \'' + var.long_name + '\' at time \'' + str(time) + '\'')

            z_index = 0
            for depth in self.points.read_axis_z():
                data = self.points.read_variable_sea_water_turbidity_at_time_and_depth(time, depth)

                var[time_index:time_index + 1, z_index:z_index + 1] = data
                z_index = z_index + 1

            time_index += 1

    def write_variable_sea_water_electrical_conductivity(self):
        var = self.ncfile.createVariable('sea_water_electrical_conductivity', float32, ('time','depth', 'point',),
                                         fill_value=9.96921e+36)
        var.long_name = "Sea Water Electrical Conductivity";
        var.standard_name = "sea_water_electrical_conductivity";
        var.units = VariableUnits.CANONICAL_UNITS[var.standard_name];

        logging.info('[DefaultTimeLevelMultiPointWriter] Writing variable \'' + var.long_name + '\'')

        time_index = 0
        for time in self.points.read_axis_t():
            logging.debug(
                '[DefaultTimeLevelMultiPointWriter] Writing variable \'' + var.long_name + '\' at time \'' + str(time) + '\'')

            z_index = 0
            for depth in self.points.read_axis_z():
                data = self.points.read_variable_sea_water_conductivity_at_time_and_depth(time, depth)

                var[time_index:time_index + 1, z_index:z_index + 1] = data
                z_index = z_index + 1

            time_index += 1