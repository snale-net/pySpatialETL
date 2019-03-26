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
from utils.VariableDefinition import VariableDefinition
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

        # Time dimension
        self.ncfile.createDimension(VariableDefinition.VARIABLE_NAME['time'], None)
        times = self.ncfile.createVariable('time', float64, (VariableDefinition.VARIABLE_NAME['time'],))
        times.units= 'seconds since 1970-01-01 00:00:00'
        times.calendar= 'gregorian'
        times.standard_name= 'time'
        times.axis='T'
        times.conventions = "UTC time"

        times[:] = date2num(self.points.read_axis_t(), units = times.units, calendar = times.calendar)

        # Depth dimension
        self.ncfile.createDimension(VariableDefinition.VARIABLE_NAME['depth'], self.points.get_z_size())
        var = self.ncfile.createVariable(VariableDefinition.VARIABLE_NAME['depth'], float64,
                                         (VariableDefinition.VARIABLE_NAME['depth'],))
        var.standard_name = VariableDefinition.STANDARD_NAME['depth']
        var.long_name = VariableDefinition.LONG_NAME['depth']
        var.positive = "down";
        var.axis = 'Z'
        var.units = VariableDefinition.CANONICAL_UNITS['depth'];
        var[:] = self.points.read_axis_z()


    def close(self):
        self.ncfile.close()

    def write_variable_baroclinic_sea_water_velocity(self):

        ucur = self.ncfile.createVariable(VariableDefinition.VARIABLE_NAME['baroclinic_eastward_sea_water_velocity'], float32, (VariableDefinition.VARIABLE_NAME['time'],VariableDefinition.VARIABLE_NAME['depth'],VariableDefinition.VARIABLE_NAME['point'],),
                                          fill_value=9.96921e+36)
        ucur.long_name = VariableDefinition.LONG_NAME['baroclinic_eastward_sea_water_velocity']
        ucur.standard_name = VariableDefinition.STANDARD_NAME['baroclinic_eastward_sea_water_velocity']
        ucur.units = VariableDefinition.CANONICAL_UNITS['baroclinic_eastward_sea_water_velocity']
        ucur.comment = "cur=sqrt(U**2+V**2)";

        vcur = self.ncfile.createVariable(VariableDefinition.VARIABLE_NAME['baroclinic_northward_sea_water_velocity'], float32, (VariableDefinition.VARIABLE_NAME['time'],VariableDefinition.VARIABLE_NAME['depth'],VariableDefinition.VARIABLE_NAME['point']),
                                          fill_value=9.96921e+36)
        vcur.long_name = VariableDefinition.LONG_NAME['baroclinic_northward_sea_water_velocity']
        vcur.standard_name = VariableDefinition.STANDARD_NAME['baroclinic_northward_sea_water_velocity']
        vcur.units = VariableDefinition.CANONICAL_UNITS['baroclinic_northward_sea_water_velocity']
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

        var = self.ncfile.createVariable(VariableDefinition.VARIABLE_NAME['sea_water_temperature'], float32, (VariableDefinition.VARIABLE_NAME['time'],VariableDefinition.VARIABLE_NAME['depth'],VariableDefinition.VARIABLE_NAME['point']), fill_value=9.96921e+36)
        var.long_name = VariableDefinition.LONG_NAME['sea_water_temperature']
        var.standard_name = VariableDefinition.STANDARD_NAME['sea_water_temperature']
        var.units = VariableDefinition.CANONICAL_UNITS['sea_water_temperature']

        logging.info('[DefaultTimeLevelMultiPointWriter] Writing variable \'' + str(VariableDefinition.LONG_NAME['sea_water_temperature']) + '\'')

        time_index = 0
        for time in self.points.read_axis_t():
            logging.debug('[DefaultTimeLevelMultiPointWriter] Writing variable \''+str(VariableDefinition.LONG_NAME['sea_water_temperature'])+'\' at time \'' + str(time) + '\'')

            z_index = 0
            for depth in self.points.read_axis_z():
                data = self.points.read_variable_sea_water_temperature_at_time_and_depth(time,depth)

                var[time_index:time_index + 1,z_index:z_index + 1] = data
                z_index = z_index +1

            time_index += 1

    def write_variable_sea_water_salinity(self):
        var = self.ncfile.createVariable(VariableDefinition.VARIABLE_NAME['sea_water_salinity'], float32, (VariableDefinition.VARIABLE_NAME['time'],VariableDefinition.VARIABLE_NAME['depth'],VariableDefinition.VARIABLE_NAME['point']), fill_value=9.96921e+36)
        var.long_name = VariableDefinition.LONG_NAME['sea_water_salinity']
        var.standard_name = VariableDefinition.STANDARD_NAME['sea_water_salinity']
        var.units = VariableDefinition.CANONICAL_UNITS['sea_water_salinity']

        logging.info('[DefaultTimeLevelMultiPointWriter] Writing variable \'' + str(VariableDefinition.LONG_NAME['sea_water_salinity']) + '\'')

        time_index = 0
        for time in self.points.read_axis_t():
            logging.debug(
                '[DefaultTimeLevelMultiPointWriter] Writing variable \'' +str(VariableDefinition.LONG_NAME['sea_water_salinity']) + '\' at time \'' + str(time) + '\'')

            z_index = 0
            for depth in self.points.read_axis_z():
                data = self.points.read_variable_sea_water_salinity_at_time_and_depth(time, depth)

                var[time_index:time_index + 1, z_index:z_index + 1] = data
                z_index = z_index + 1

            time_index += 1

    def write_variable_sea_water_density(self):
        var = self.ncfile.createVariable(VariableDefinition.VARIABLE_NAME['sea_water_density'], float32, (VariableDefinition.VARIABLE_NAME['time'],VariableDefinition.VARIABLE_NAME['depth'],VariableDefinition.VARIABLE_NAME['point']), fill_value=9.96921e+36)
        var.long_name = VariableDefinition.LONG_NAME['sea_water_density']
        var.standard_name = VariableDefinition.STANDARD_NAME['sea_water_density']
        var.units = VariableDefinition.CANONICAL_UNITS['sea_water_density']

        logging.info('[DefaultTimeLevelMultiPointWriter] Writing variable \'' + str(VariableDefinition.LONG_NAME['sea_water_density']) + '\'')

        time_index = 0
        for time in self.points.read_axis_t():
            logging.debug(
                '[DefaultTimeLevelMultiPointWriter] Writing variable \'' + str(VariableDefinition.LONG_NAME['sea_water_density']) + '\' at time \'' + str(time) + '\'')

            z_index = 0
            for depth in self.points.read_axis_z():
                data = self.points.read_variable_sea_water_density_at_time_and_depth(time, depth)

                var[time_index:time_index + 1, z_index:z_index + 1] = data
                z_index = z_index + 1

            time_index += 1

    def write_variable_sea_water_turbidity(self):
        var = self.ncfile.createVariable(VariableDefinition.VARIABLE_NAME['sea_water_turbidity'], float32, (VariableDefinition.VARIABLE_NAME['time'],VariableDefinition.VARIABLE_NAME['depth'],VariableDefinition.VARIABLE_NAME['point']), fill_value=9.96921e+36)
        var.long_name = VariableDefinition.LONG_NAME['sea_water_turbidity']
        var.standard_name = VariableDefinition.STANDARD_NAME['sea_water_turbidity']
        var.units = VariableDefinition.CANONICAL_UNITS['sea_water_turbidity']

        logging.info('[DefaultTimeLevelMultiPointWriter] Writing variable \'' + str(VariableDefinition.LONG_NAME['sea_water_turbidity']) + '\'')

        time_index = 0
        for time in self.points.read_axis_t():
            logging.debug(
                '[DefaultTimeLevelMultiPointWriter] Writing variable \'' +  str(VariableDefinition.LONG_NAME['sea_water_turbidity']) + '\' at time \'' + str(time) + '\'')

            z_index = 0
            for depth in self.points.read_axis_z():
                data = self.points.read_variable_sea_water_turbidity_at_time_and_depth(time, depth)

                var[time_index:time_index + 1, z_index:z_index + 1] = data
                z_index = z_index + 1

            time_index += 1

    def write_variable_sea_water_electrical_conductivity(self):
        var = self.ncfile.createVariable(VariableDefinition.VARIABLE_NAME['sea_water_electrical_conductivity'], float32, (VariableDefinition.VARIABLE_NAME['time'],VariableDefinition.VARIABLE_NAME['depth'],VariableDefinition.VARIABLE_NAME['point']),
                                         fill_value=9.96921e+36)
        var.long_name = VariableDefinition.LONG_NAME['sea_water_electrical_conductivity']
        var.standard_name = VariableDefinition.STANDARD_NAME['sea_water_electrical_conductivity']
        var.units = VariableDefinition.CANONICAL_UNITS['sea_water_electrical_conductivity']

        logging.info('[DefaultTimeLevelMultiPointWriter] Writing variable \'' + str(VariableDefinition.LONG_NAME['sea_water_electrical_conductivity']) + '\'')

        time_index = 0
        for time in self.points.read_axis_t():
            logging.debug(
                '[DefaultTimeLevelMultiPointWriter] Writing variable \'' + str(VariableDefinition.LONG_NAME['sea_water_electrical_conductivity']) + '\' at time \'' + str(time) + '\'')

            z_index = 0
            for depth in self.points.read_axis_z():
                data = self.points.read_variable_sea_water_conductivity_at_time_and_depth(time, depth)

                var[time_index:time_index + 1, z_index:z_index + 1] = data
                z_index = z_index + 1

            time_index += 1