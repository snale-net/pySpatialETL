# #! /usr/bin/env python2.7
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
from spatialetl.coverage.io.CoverageWriter import CoverageWriter
from netCDF4 import Dataset
from netCDF4 import date2num
from numpy import float32
from numpy import float64
import numpy as np
import logging

class WW3Writer (CoverageWriter):

    def __init__(self, cov,myFile,depths):
        CoverageWriter.__init__(self,cov,myFile);
        self.ncfile = None

        if self.coverage.is_regular_grid()==True:
            raise IOError("This writer is specific to non-regular grid.")   
        
        self.ncfile = Dataset(self.filename, 'w', format='NETCDF4')
        self.ncfile.description = 'Generated with pyGeoSpatialETL'

        # dimensions
        self.ncfile.createDimension('time', None)
        self.ncfile.createDimension('latitude', self.coverage.get_y_size())
        self.ncfile.createDimension('longitude', self.coverage.get_x_size())

        # variables
        times = self.ncfile.createVariable('time', float64, ('time',))
        times.units= 'seconds since 1970-01-01 00:00:00' 
        times.calendar= 'gregorian'
        times.standard_name= 'time'
        times.axis='T'
        times.conventions = "UTC time"
        
        latitudes = self.ncfile.createVariable('latitude', float32, ('latitude','longitude',))
        latitudes.units = "degree_north" ;
        latitudes.long_name = "latitude" ;
        latitudes.standard_name = "latitude" ;
        latitudes.valid_min = "-90.0";
        latitudes.valid_max = "90.0" ;
        latitudes.axis = "Y" ;
        
        longitudes = self.ncfile.createVariable('longitude', float32, ('latitude','longitude',))
        longitudes.units = "degree_east" ;
        longitudes.long_name = "longitude" ;
        longitudes.standard_name = "longitude" ;
        longitudes.valid_min = "-180.0" ;
        longitudes.valid_max = "180.0" ;
        longitudes.axis = "X" ; 
        
         # data
        latitudes[:,:] = self.coverage.read_axis_y();
        longitudes[:,:] = self.coverage.read_axis_x();
        times[:] = date2num(self.coverage.read_axis_t(), units = times.units, calendar = times.calendar)

    def close(self):
        self.ncfile.close()
            
    def write_variable_ssh(self):
        if self.ncfile == None:
            raise IOError("Please call write_axis() first") 
            
        wlv = self.ncfile.createVariable('wlv', float32, ('time', 'latitude', 'longitude',),fill_value=9.96921e+36)
        wlv.long_name = "sea surface height above sea level" ;
        wlv.standard_name = "sea_surface_height_above_sea_level" ;
        wlv.globwave_name = "sea_surface_height_above_sea_level" ;
        wlv.units = "m" ;
        
        time_index=0
        for time in self.coverage.read_axis_t():
            logging.info('[WW3Writer] Writing variable \'wlv\' at time \''+str(time)+'\'')
            wlv[time_index:time_index+1,:] = self.coverage.read_variable_ssh_at_time(time)
            #s = self.coverage.read_variable_ssh_at_time(time)
            #s += 0.4466
            #wlv[time_index:time_index+1,:] = s
            time_index += 1
            
    def write_variable_current_at_depth(self,z):
        
        if self.ncfile == None:
            raise IOError("Please call write_axis() first")   
            
        ucur = self.ncfile.createVariable('ucur', float32, ('time', 'latitude', 'longitude',),fill_value=9.96921e+36)
        ucur.long_name = "eastward current" ;
        ucur.standard_name = "eastward_sea_water_velocity" ;
        ucur.globwave_name = "eastward_sea_water_velocity" ;
        ucur.units = "m s-1" ;
        #ucur.scale_factor = 1.f ;
        #ucur.add_offset = 0.f ;
        #ucur.valid_min = -990 ;
        #ucur.valid_max = 990 ;
        ucur.comment = "cur=sqrt(U**2+V**2)" ;
        
        vcur = self.ncfile.createVariable('vcur', float32, ('time', 'latitude', 'longitude',),fill_value=9.96921e+36)
        vcur.long_name = "northward current" ;
        vcur.standard_name = "northward_sea_water_velocity" ;
        vcur.globwave_name = "northward_sea_water_velocity" ;
        vcur.units = "m s-1" ;
        #ucur.scale_factor = 1.f ;
        #ucur.add_offset = 0.f ;
        #ucur.valid_min = -990 ;
        #ucur.valid_max = 990 ;
        vcur.comment = "cur=sqrt(U**2+V**2)" ;

        time_index=0
        for time in self.coverage.read_axis_t():

            logging.info('[WW3Writer] Writing variable \'current\' at time \''+str(time)+'\'')
            cur = self.coverage.read_variable_current_at_time_and_depth(time,z)

            ucur[time_index:time_index+1,:,:] = cur[0]
            vcur[time_index:time_index+1,:,:] = cur[1]
            time_index += 1

    def write_variable_current(self):

        if self.ncfile == None:
            raise IOError("Please call write_axis() first")

        ucur = self.ncfile.createVariable('ucur', float32, ('time', 'latitude', 'longitude',),fill_value=9.96921e+36)
        ucur.long_name = "eastward current" ;
        ucur.standard_name = "eastward_sea_water_velocity" ;
        ucur.globwave_name = "eastward_sea_water_velocity" ;
        ucur.units = "m s-1" ;
        #ucur.scale_factor = 1.f ;
        #ucur.add_offset = 0.f ;
        #ucur.valid_min = -990 ;
        #ucur.valid_max = 990 ;
        ucur.comment = "cur=sqrt(U**2+V**2)" ;

        vcur = self.ncfile.createVariable('vcur', float32, ('time', 'latitude', 'longitude',),fill_value=9.96921e+36)
        vcur.long_name = "northward current" ;
        vcur.standard_name = "northward_sea_water_velocity" ;
        vcur.globwave_name = "northward_sea_water_velocity" ;
        vcur.units = "m s-1" ;
        #ucur.scale_factor = 1.f ;
        #ucur.add_offset = 0.f ;
        #ucur.valid_min = -990 ;
        #ucur.valid_max = 990 ;
        vcur.comment = "cur=sqrt(U**2+V**2)" ;

        time_index=0
        for time in self.coverage.read_axis_t():

            logging.info('[WW3Writer] Writing variable \'current\' at time \''+str(time)+'\'')
            cur = self.coverage.read_variable_current_at_time(time)

            ucur[time_index:time_index+1,:,:] = cur[0]
            vcur[time_index:time_index+1,:,:] = cur[1]
            time_index += 1


    def write_variable_hs(self):
        if self.ncfile == None:
            raise IOError("Please call write_axis() first")

        var = self.ncfile.createVariable('hs', float32, ('time', 'latitude', 'longitude',),fill_value=9.96921e+36)
        var.long_name = "sea surface wave height" ;
        var.standard_name = "sea_surface_wave_height" ;
        var.globwave_name = "sea_surface_wave_height" ;
        var.units = "m" ;
        #wlv.scale_factor = "1.f" ;
        #wlv.add_offset = "0.f" ;
        #wlv.valid_min = "0f" ;
        #wlv.valid_max = 10000f ;

        time_index=0
        for time in self.coverage.read_axis_t():
            logging.info('[WW3Writer] Writing variable \'hs\' at time \''+str(time)+'\'')
            var[time_index:time_index+1,:] = self.coverage.read_variable_hs_at_time(time)
            time_index += 1
        
    def write_variable_2D_mask(self):
        if self.ncfile == None:
            raise IOError("Please call write_axis() first")

        var = self.ncfile.createVariable('mask', int, ( 'latitude', 'longitude',),fill_value=-999)
        var.long_name = "land/sea mask" ;
        var.standard_name = "land_sea_mask" ;
        var.globwave_name = "land_sea_mask" ;
        var.units = "m" ;
        #wlv.scale_factor = "1.f" ;
        #wlv.add_offset = "0.f" ;
        #wlv.valid_min = "0f" ;
        #wlv.valid_max = 10000f ;

        logging.info('[WW3Writer] Writing variable \'2D mask\'')
        var[:] = self.coverage.read_variable_2D_mask()

    def write_variable_waves_dir(self):
        if self.ncfile == None:
            raise IOError("Please call write_axis() first")

        var = self.ncfile.createVariable('waves_dir', float32, ('time', 'latitude', 'longitude',),fill_value=9.96921e+36)
        var.long_name = "waves_dir" ;
        var.standard_name = "waves_dir" ;
        var.globwave_name = "waves_dir" ;
        var.units = "deg" ;
        #wlv.scale_factor = "1.f" ;
        #wlv.add_offset = "0.f" ;
        #wlv.valid_min = "0f" ;
        #wlv.valid_max = 10000f ;

        time_index=0
        for time in self.coverage.read_axis_t():
            logging.info('[WW3Writer] Writing variable \'waves_dir\' at time \''+str(time)+'\'')
            var[time_index:time_index+1,:] = self.coverage.read_variable_waves_dir_at_time(time)
            time_index += 1

    def write_variable_surface_stokes_drift(self):
        if self.ncfile == None:
            raise IOError("Please call write_axis() first")

        ucur = self.ncfile.createVariable('uss', float32, ('time', 'latitude', 'longitude',),fill_value=9.96921e+36)
        ucur.long_name = "eastward surface_stokes_drift" ;
        ucur.standard_name = "eastward_sea_water_velocity" ;
        ucur.globwave_name = "eastward_sea_water_velocity" ;
        ucur.units = "m s-1" ;
        #ucur.scale_factor = 1.f ;
        #ucur.add_offset = 0.f ;
        #ucur.valid_min = -990 ;
        #ucur.valid_max = 990 ;
        ucur.comment = "cur=sqrt(U**2+V**2)" ;

        vcur = self.ncfile.createVariable('vss', float32, ('time', 'latitude', 'longitude',),fill_value=9.96921e+36)
        vcur.long_name = "northward surface_stokes_drift" ;
        vcur.standard_name = "northward_sea_water_velocity" ;
        vcur.globwave_name = "northward_sea_water_velocity" ;
        vcur.units = "m s-1" ;
        #ucur.scale_factor = 1.f ;
        #ucur.add_offset = 0.f ;
        #ucur.valid_min = -990 ;
        #ucur.valid_max = 990 ;
        vcur.comment = "cur=sqrt(U**2+V**2)" ;

        time_index=0
        for time in self.coverage.read_axis_t():

            logging.info('[WW3Writer] Writing variable \'surface_stokes_drift\' at time \''+str(time)+'\'')
            cur = self.coverage.read_variable_surface_stokes_drift_at_time(time)

            ucur[time_index:time_index+1,:,:] = cur[0]
            vcur[time_index:time_index+1,:,:] = cur[1]
            time_index += 1

       
    