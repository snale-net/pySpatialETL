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

from coverage.operator.interpolator.InterpolatorCore import resample_2d_to_grid
from coverage import TimeCoverage
from coverage.LevelCoverage import LevelCoverage
from coverage.TimeLevelCoverage import TimeLevelCoverage
import logging
from coverage.io.CoverageWriter import CoverageWriter
from netCDF4 import Dataset
from netCDF4 import date2num
from numpy import float32
from numpy import float64
import numpy as np

class CoverageInterpolator(CoverageWriter):
    """
    Cette classe permet d'interpoler un coverage non régulier sur une grille régulière.
    Cette classe écrit un fichier netcdf
    """

    def __init__(self, cov,resX,resY,myFile, depths):
        """
    Constructeur
    @param cov : la coverage
    @param resX : résolution souhaitée en X, en degrès
    @param resY : résolution souhaitée en Y, en degrès
    @param myFile : fichier de destination
    @param depths : tableau des profondeurs souhaitée (entier = couche, flottant = profondeur)
    """
        CoverageWriter.__init__(self,cov,myFile);

        self.targetResX = resX
        self.targetResY = resY;
        self.targetDepths = depths

        # we compute the destination grid
        Ymin=np.min(self.coverage.read_axis_y())
        Ymax=np.max(self.coverage.read_axis_y())
        Xmin=np.min(self.coverage.read_axis_x())
        Xmax=np.max(self.coverage.read_axis_x())

        res=np.mean([self.targetResX,self.targetResY])
        self.lon_reg,self.lat_reg=np.meshgrid(np.arange(Xmin, Xmax, res),np.arange(Ymin, Ymax, res))
        targetAxisX = self.lon_reg[0,:]
        targetAxisY= self.lat_reg[:,0]

        self.ncfile = Dataset(self.filename, 'w', format='NETCDF4')
        self.ncfile.description = 'Coverage Interpolator. Generated with Coverage Processing tools'

        # dimensions
        self.ncfile.createDimension('latitude', np.shape(targetAxisY)[0])
        self.ncfile.createDimension('longitude', np.shape(targetAxisX)[0])

        # variables
        latitudes = self.ncfile.createVariable('latitude', float32, ('latitude',))
        latitudes.units = "degree_north" ;
        latitudes.long_name = "latitude" ;
        latitudes.standard_name = "latitude" ;
        latitudes.valid_min = "-90.";
        latitudes.valid_max = "90." ;
        latitudes.axis = "Y" ;

        longitudes = self.ncfile.createVariable('longitude', float32, ('longitude',))
        longitudes.units = "degree_east" ;
        longitudes.long_name = "longitude" ;
        longitudes.standard_name = "longitude" ;
        longitudes.valid_min = "-180." ;
        longitudes.valid_max = "180." ;
        longitudes.axis = "X" ;

        # data
        latitudes[:] = targetAxisY;
        longitudes[:] = targetAxisX;

        if(isinstance(self.coverage, TimeCoverage) or isinstance(self.coverage, TimeLevelCoverage)):

            self.ncfile.createDimension('time', None)
            times = self.ncfile.createVariable('time', float64, ('time',))
            times.units= 'seconds since 1970-01-01 00:00:00'
            times.calendar= 'gregorian'
            times.standard_name= 'time'
            times.axis='T'
            times.conventions = "UTC time"

            times[:] = date2num(self.coverage.read_axis_t(), units = times.units, calendar = times.calendar)

        if(isinstance(self.coverage, LevelCoverage) or isinstance(self.coverage, TimeLevelCoverage)):

            self.ncfile.createDimension('depth', np.size(self.targetDepths))
            levels = self.ncfile.createVariable('depth', float64, ('depth',))
            levels.standard_name= 'depth'
            levels.long_name="Positive depth"
            levels.units = "m" ;
            levels.axis='Z'

            levels[:] = self.targetDepths

    def close(self):
        self.ncfile.close()

    def resample_variable_topography(self):
        if self.ncfile == None:
            raise IOError("Please call write_axis() first")

        logging.info('[CoverageInterpolator] Resample variable \'topography\' to the resolution '+str(self.targetResX)+'/'+str(self.targetResY)+'.')

        var = self.ncfile.createVariable('topo', float32, ('latitude', 'longitude',),fill_value=9.96921e+36)
        var.long_name = "topography" ;
        var.standard_name = "topography" ;
        var.globwave_name = "topography" ;
        var.units = "m" ;
        #wlv.scale_factor = "1.f" ;
        #wlv.add_offset = "0.f" ;
        #wlv.valid_min = "0f" ;
        #wlv.valid_max = 10000f ;

        var[:,:] = resample_2d_to_grid(self.coverage.read_axis_x(),self.coverage.read_axis_y(),self.lon_reg,self.lat_reg,self.coverage.read_variable_topography())

    def resample_variable_bathymetry(self):
        if self.ncfile == None:
            raise IOError("Please call write_axis() first")

        logging.info('[CoverageInterpolator] Resample variable \'bathymetry\' to the resolution '+str(self.targetResX)+'/'+str(self.targetResY)+'.')

        var = self.ncfile.createVariable('bathy', float32, ('latitude', 'longitude',),fill_value=9.96921e+36)
        var.long_name = "bathymetry" ;
        var.standard_name = "bathymetry" ;
        var.globwave_name = "bathymetry" ;
        var.units = "m" ;
        #wlv.scale_factor = "1.f" ;
        #wlv.add_offset = "0.f" ;
        #wlv.valid_min = "0f" ;
        #wlv.valid_max = 10000f ;

        var[:,:] = resample_2d_to_grid(self.coverage.read_axis_x(),self.coverage.read_axis_y(),self.lon_reg,self.lat_reg,self.coverage.read_variable_bathymetry())

    def resample_variable_mesh_size(self):
        if self.ncfile == None:
            raise IOError("Please call write_axis() first")

        logging.info('[CoverageInterpolator] Resample variable \'mesh_size\' to the resolution '+str(self.targetResX)+'/'+str(self.targetResY)+'.')

        var = self.ncfile.createVariable('mesh_size', float32, ('latitude', 'longitude',),fill_value=9.96921e+36)
        var.long_name = "mesh_size" ;
        var.standard_name = "mesh_size" ;
        var.globwave_name = "mesh_size" ;
        var.units = "m" ;
        #wlv.scale_factor = "1.f" ;
        #wlv.add_offset = "0.f" ;
        #wlv.valid_min = "0f" ;
        #wlv.valid_max = 10000f ;

        var[:,:] = resample_2d_to_grid(self.coverage.read_axis_x(),self.coverage.read_axis_y(),self.lon_reg,self.lat_reg,self.coverage.read_variable_mesh_size())

    def resample_variable_Ha(self):
        if self.ncfile == None:
            raise IOError("Please call write_axis() first")

        logging.info('[CoverageInterpolator] Resample variable \'Ha\' to the resolution '+str(self.targetResX)+'/'+str(self.targetResY)+'.')

        var = self.ncfile.createVariable('Ha', float32, ('latitude', 'longitude',),fill_value=9.96921e+36)
        var.long_name = "Amplitude tide comp" ;
        var.standard_name = "amp_tide_comp" ;
        var.globwave_name = "amp_tide_comp" ;
        var.units = "m" ;
        #wlv.scale_factor = "1.f" ;
        #wlv.add_offset = "0.f" ;
        #wlv.valid_min = "0f" ;
        #wlv.valid_max = 10000f ;

        var[:,:] = resample_2d_to_grid(self.coverage.read_axis_x(),self.coverage.read_axis_y(),self.lon_reg,self.lat_reg,self.coverage.read_variable_Ha())

    def resample_variable_waves_bottom_dissipation(self):
        if self.ncfile == None:
            raise IOError("Please call write_axis() first")

        logging.info('[CoverageInterpolator] Resample variable \'waves_bottom_dissipation\' at resolution '+str(self.targetResX)+'/'+str(self.targetResY)+'.')

        wlv = self.ncfile.createVariable('waves_bottom_dissipation', float32, ('time', 'latitude', 'longitude',),fill_value=9.96921e+36)
        wlv.long_name = "waves bottom dissipation" ;
        wlv.standard_name = "waves_bottom_dissipation" ;
        wlv.globwave_name = "waves_bottom_dissipation" ;
        wlv.units = "W m-2" ;
        #wlv.scale_factor = "1.f" ;
        #wlv.add_offset = "0.f" ;
        #wlv.valid_min = "0f" ;
        #wlv.valid_max = 10000f ;

        time_index=0
        for time in self.coverage.read_axis_t():
            wlv[time_index:time_index+1,:,:] = resample_2d_to_grid(self.coverage.read_axis_x(),self.coverage.read_axis_y(),self.lon_reg,self.lat_reg,self.coverage.read_variable_waves_bottom_dissipation_at_time(time))
            time_index += 1

    def resample_variable_ssh(self):
        if self.ncfile == None:
            raise IOError("Please call write_axis() first")

        logging.info('[CoverageInterpolator] Resample variable \'ssh\' at resolution '+str(self.targetResX)+'/'+str(self.targetResY)+'.')

        wlv = self.ncfile.createVariable('ssh', float32, ('time', 'latitude', 'longitude',),fill_value=9.96921e+36)
        wlv.long_name = "sea surface height above sea level" ;
        wlv.standard_name = "sea_surface_height_above_sea_level" ;
        wlv.globwave_name = "sea_surface_height_above_sea_level" ;
        wlv.units = "m" ;
        #wlv.scale_factor = "1.f" ;
        #wlv.add_offset = "0.f" ;
        #wlv.valid_min = "0f" ;
        #wlv.valid_max = 10000f ;

        time_index=0
        for time in self.coverage.read_axis_t():
            wlv[time_index:time_index+1,:,:] = resample_2d_to_grid(self.coverage.read_axis_x(),self.coverage.read_axis_y(),self.lon_reg,self.lat_reg,self.coverage.read_variable_ssh_at_time(time))
            time_index += 1

    def resample_variable_bathy_ssh(self):
        if self.ncfile == None:
            raise IOError("Please call write_axis() first")

        logging.info('[CoverageInterpolator] Resample variable \'bathy_ssh\' at resolution '+str(self.targetResX)+'/'+str(self.targetResY)+'.')

        wlv = self.ncfile.createVariable('bathy_ssh', float32, ('time', 'latitude', 'longitude',),fill_value=9.96921e+36)
        wlv.long_name = "sea surface height above sea level" ;
        wlv.standard_name = "sea_surface_height_above_sea_level" ;
        wlv.globwave_name = "sea_surface_height_above_sea_level" ;
        wlv.units = "m" ;
        #wlv.scale_factor = "1.f" ;
        #wlv.add_offset = "0.f" ;
        #wlv.valid_min = "0f" ;
        #wlv.valid_max = 10000f ;

        time_index=0
        for time in self.coverage.read_axis_t():
            wlv[time_index:time_index+1,:,:] = resample_2d_to_grid(self.coverage.read_axis_x(),self.coverage.read_axis_y(),self.lon_reg,self.lat_reg,self.coverage.read_variable_bathy_ssh_at_time(time))
            time_index += 1

    def resample_variable_baroclinic_sea_water_velocity(self,vertical_method="nearest"):
        """
        Interpole un champ de courant au niveau donnée à la construction
        @vertical_method: méthode d'interpolation verticale. "nearest" ou "linear"
        """
        if self.ncfile == None:
            raise IOError("Please call write_axis() first")

        ucur = self.ncfile.createVariable('baroclinic_eastward_sea_water_velocity', float32, ('time', 'depth', 'latitude', 'longitude',),fill_value=9.96921e+36)
        ucur.long_name = "baroclinic_eastward_sea_water_velocity" ;
        ucur.standard_name = "eastward_sea_water_velocity" ;
        ucur.globwave_name = "eastward_sea_water_velocity" ;
        ucur.units = "m s-1" ;
        #ucur.scale_factor = 1.f ;
        #ucur.add_offset = 0.f ;
        #ucur.valid_min = -990 ;
        #ucur.valid_max = 990 ;
        ucur.comment = "cur=sqrt(U**2+V**2)" ;

        vcur = self.ncfile.createVariable('baroclinic_northward_sea_water_velocity', float32, ('time', 'depth', 'latitude', 'longitude',),fill_value=9.96921e+36)
        vcur.long_name = "baroclinic_northward_sea_water_velocity" ;
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

            logging.info('[CoverageInterpolator] Resample variable \'Baroclinic Current\' at time '+str(time)+' to the resolution '+str(self.targetResX)+'/'+str(self.targetResY)+'.')

            level_index = 0
            for level in self.targetDepths:

                if type(level) == int:
                    logging.info('[CoverageInterpolator] At index level '+str(level)+'')
                else:
                    logging.info('[CoverageInterpolator] At depth '+str(level)+' m.')

                cur = self.coverage.read_variable_current_at_time_and_depth(time,level)
                ucur[time_index:time_index+1,level_index:level_index+1,:,:] = resample_2d_to_grid(self.coverage.read_axis_x(),self.coverage.read_axis_y(),self.lon_reg,self.lat_reg,cur[0])
                vcur[time_index:time_index+1,level_index:level_index+1,:,:] = resample_2d_to_grid(self.coverage.read_axis_x(),self.coverage.read_axis_y(),self.lon_reg,self.lat_reg,cur[1])
                level_index+= 1

            time_index += 1

    def resample_variable_barotropic_sea_water_velocity(self):

        if self.ncfile == None:
            raise IOError("Please call write_axis() first")

        logging.info('[CoverageInterpolator] Resample variable \'Barotropic Current\' to the resolution '+str(self.targetResX)+'/'+str(self.targetResY)+'.')

        ucur = self.ncfile.createVariable('barotropic_eastward_sea_water_velocity', float32, ('time', 'latitude', 'longitude',),fill_value=9.96921e+36)
        ucur.long_name = "barotropic_eastward_sea_water_velocity" ;
        ucur.standard_name = "eastward_sea_water_velocity" ;
        ucur.units = "m s-1" ;
        #ucur.scale_factor = 1.f ;
        #ucur.add_offset = 0.f ;
        #ucur.valid_min = -990 ;
        #ucur.valid_max = 990 ;
        ucur.comment = "cur=sqrt(U**2+V**2)" ;

        vcur = self.ncfile.createVariable('barotropic_northward_sea_water_velocity', float32, ('time', 'latitude', 'longitude',),fill_value=9.96921e+36)
        vcur.long_name = "barotropic_northward_sea_water_velocity" ;
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
            cur = self.coverage.read_variable_current_at_time(time)
            ucur[time_index:time_index+1,:,:] = resample_2d_to_grid(self.coverage.read_axis_x(),self.coverage.read_axis_y(),self.lon_reg,self.lat_reg,cur[0])
            vcur[time_index:time_index+1,:,:] = resample_2d_to_grid(self.coverage.read_axis_x(),self.coverage.read_axis_y(),self.lon_reg,self.lat_reg,cur[1])
            time_index += 1

    def resample_variable_taw(self):

        if self.ncfile == None:
            raise IOError("Please call write_axis() first")

        logging.info('[CoverageInterpolator] Resample variable \'taw\' to the resolution '+str(self.targetResX)+'/'+str(self.targetResY)+'.')

        ucur = self.ncfile.createVariable('utaw', float32, ('time', 'latitude', 'longitude',),fill_value=9.96921e+36)
        ucur.long_name = "eastward current" ;
        ucur.standard_name = "eastward_sea_water_velocity" ;
        ucur.globwave_name = "eastward_sea_water_velocity" ;
        ucur.units = "m s-1" ;
        #ucur.scale_factor = 1.f ;
        #ucur.add_offset = 0.f ;
        #ucur.valid_min = -990 ;
        #ucur.valid_max = 990 ;
        ucur.comment = "cur=sqrt(U**2+V**2)" ;

        vcur = self.ncfile.createVariable('vtaw', float32, ('time', 'latitude', 'longitude',),fill_value=9.96921e+36)
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
            cur = self.coverage.read_variable_taw_at_time(time)
            ucur[time_index:time_index+1,:,:] = resample_2d_to_grid(self.coverage.read_axis_x(),self.coverage.read_axis_y(),self.lon_reg,self.lat_reg,cur[0])
            vcur[time_index:time_index+1,:,:] = resample_2d_to_grid(self.coverage.read_axis_x(),self.coverage.read_axis_y(),self.lon_reg,self.lat_reg,cur[1])
            time_index += 1

    def resample_variable_two(self):

        if self.ncfile == None:
            raise IOError("Please call write_axis() first")

        logging.info('[CoverageInterpolator] Resample variable \'two\' to the resolution '+str(self.targetResX)+'/'+str(self.targetResY)+'.')

        ucur = self.ncfile.createVariable('utwo', float32, ('time', 'latitude', 'longitude',),fill_value=9.96921e+36)
        ucur.long_name = "eastward current" ;
        ucur.standard_name = "eastward_sea_water_velocity" ;
        ucur.globwave_name = "eastward_sea_water_velocity" ;
        ucur.units = "m s-1" ;
        #ucur.scale_factor = 1.f ;
        #ucur.add_offset = 0.f ;
        #ucur.valid_min = -990 ;
        #ucur.valid_max = 990 ;
        ucur.comment = "cur=sqrt(U**2+V**2)" ;

        vcur = self.ncfile.createVariable('vtwo', float32, ('time', 'latitude', 'longitude',),fill_value=9.96921e+36)
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
            cur = self.coverage.read_variable_two_at_time(time)
            ucur[time_index:time_index+1,:,:] = resample_2d_to_grid(self.coverage.read_axis_x(),self.coverage.read_axis_y(),self.lon_reg,self.lat_reg,cur[0])
            vcur[time_index:time_index+1,:,:] = resample_2d_to_grid(self.coverage.read_axis_x(),self.coverage.read_axis_y(),self.lon_reg,self.lat_reg,cur[1])
            time_index += 1

    def resample_variable_hs(self):
        if self.ncfile == None:
            raise IOError("Please call write_axis() first")

        logging.info('[CoverageInterpolator] Resample variable \'hs\' at resolution '+str(self.targetResX)+'/'+str(self.targetResY)+'.')

        var = self.ncfile.createVariable('hs', float32, ('time', 'latitude', 'longitude',),fill_value=9.96921e+36)
        var.long_name = "surface wave height" ;
        var.standard_name = "surface_wave_height" ;
        var.globwave_name = "surface_wave_height" ;
        var.units = "m" ;
        #wlv.scale_factor = "1.f" ;
        #wlv.add_offset = "0.f" ;
        #var.valid_min = 0. ;
        #var.valid_max = 100.;

        time_index=0
        for time in self.coverage.read_axis_t():
            var[time_index:time_index+1,:,:] = resample_2d_to_grid(self.coverage.read_axis_x(),self.coverage.read_axis_y(),self.lon_reg,self.lat_reg,self.coverage.read_variable_hs_at_time(time))
            time_index += 1


    def resample_variable_waves_dir(self):
        if self.ncfile == None:
            raise IOError("Please call write_axis() first")

        logging.info('[CoverageInterpolator] Resample variable \'waves dir\' at resolution '+str(self.targetResX)+'/'+str(self.targetResY)+'.')

        var = self.ncfile.createVariable('waves_dir', float32, ('time', 'latitude', 'longitude',),fill_value=9.96921e+36)
        var.long_name = "surface wave dir" ;
        var.standard_name = "surface_wave_dir" ;
        var.globwave_name = "surface_wave_dir" ;
        var.units = "deg" ;
        #wlv.scale_factor = "1.f" ;
        #wlv.add_offset = "0.f" ;
        #wlv.valid_min = "0f" ;
        #wlv.valid_max = 10000f ;

        time_index=0
        for time in self.coverage.read_axis_t():
            var[time_index:time_index+1,:,:] = resample_2d_to_grid(self.coverage.read_axis_x(),self.coverage.read_axis_y(),self.lon_reg,self.lat_reg,self.coverage.read_variable_waves_dir_at_time(time))
            time_index += 1

    def resample_variable_waves_mean_period(self):
        if self.ncfile == None:
            raise IOError("Please call write_axis() first")

        logging.info('[CoverageInterpolator] Resample variable \'waves_mean_period\' at resolution '+str(self.targetResX)+'/'+str(self.targetResY)+'.')

        var = self.ncfile.createVariable('waves_mean_period', float32, ('time', 'latitude', 'longitude',),fill_value=9.96921e+36)
        var.long_name = "surface wave mean period" ;
        var.standard_name = "surface_wave_mean_period" ;
        var.globwave_name = "surface_wave_mean_period" ;
        var.units = "deg" ;
        #wlv.scale_factor = "1.f" ;
        #wlv.add_offset = "0.f" ;
        #wlv.valid_min = "0f" ;
        #wlv.valid_max = 10000f ;

        time_index=0
        for time in self.coverage.read_axis_t():
            var[time_index:time_index+1,:,:] = resample_2d_to_grid(self.coverage.read_axis_x(),self.coverage.read_axis_y(),self.lon_reg,self.lat_reg,self.coverage.read_variable_waves_mean_period_at_time(time))
            time_index += 1


    def resample_variable_2D_mask(self):
        if self.ncfile == None:
            raise IOError("Please call write_axis() first")

        logging.info('[CoverageInterpolator] Resample variable \'mask\' at resolution '+str(self.targetResX)+'/'+str(self.targetResY)+'.')

        var = self.ncfile.createVariable('mask', float32, ('latitude', 'longitude',),fill_value=9.96921e+36)
        var.long_name = "land/sea mask" ;
        var.standard_name = "land_sea_mask" ;
        var.globwave_name = "land_sea_mask" ;
        var.units = "m" ;
        #wlv.scale_factor = "1.f" ;
        #wlv.add_offset = "0.f" ;
        #wlv.valid_min = "0f" ;
        #wlv.valid_max = 10000f ;

        var[:,:] = resample_2d_to_grid(self.coverage.read_axis_x(),self.coverage.read_axis_y(),self.lon_reg,self.lat_reg,self.coverage.read_variable_2D_mask())

    def resample_variable_wind(self):

        if self.ncfile == None:
            raise IOError("Please call write_axis() first")

        logging.info('[CoverageInterpolator] Resample variable \'wind\' to the resolution '+str(self.targetResX)+'/'+str(self.targetResY)+'.')

        ucur = self.ncfile.createVariable('uwind', float32, ('time', 'latitude', 'longitude',),fill_value=9.96921e+36)
        ucur.long_name = "eastward wind" ;
        ucur.standard_name = "eastward_wind_velocity" ;
        ucur.globwave_name = "eastward_wind_velocity" ;
        ucur.units = "m s-1" ;
        #ucur.scale_factor = 1.f ;
        #ucur.add_offset = 0.f ;
        #ucur.valid_min = -990 ;
        #ucur.valid_max = 990 ;
        ucur.comment = "cur=sqrt(U**2+V**2)" ;

        vcur = self.ncfile.createVariable('vwind', float32, ('time', 'latitude', 'longitude',),fill_value=9.96921e+36)
        vcur.long_name = "northward wind" ;
        vcur.standard_name = "northward_wind_velocity" ;
        vcur.globwave_name = "northward_wind_velocity" ;
        vcur.units = "m s-1" ;
        #ucur.scale_factor = 1.f ;
        #ucur.add_offset = 0.f ;
        #ucur.valid_min = -990 ;
        #ucur.valid_max = 990 ;
        vcur.comment = "cur=sqrt(U**2+V**2)" ;

        time_index=0
        for time in self.coverage.read_axis_t():
            cur = self.coverage.read_variable_wind_at_time(time)
            ucur[time_index:time_index+1,:,:] = resample_2d_to_grid(self.coverage.read_axis_x(),self.coverage.read_axis_y(),self.lon_reg,self.lat_reg,cur[0])
            vcur[time_index:time_index+1,:,:] = resample_2d_to_grid(self.coverage.read_axis_x(),self.coverage.read_axis_y(),self.lon_reg,self.lat_reg,cur[1])
            time_index += 1

    def resample_variable_surface_stokes_drift(self):

        if self.ncfile == None:
            raise IOError("Please call write_axis() first")

        logging.info('[CoverageInterpolator] Resample variable \'surface_stokes_drift\' to the resolution '+str(self.targetResX)+'/'+str(self.targetResY)+'.')

        ucur = self.ncfile.createVariable('uuss', float32, ('time', 'latitude', 'longitude',),fill_value=9.96921e+36)
        ucur.long_name = "eastward surface stokes drift" ;
        ucur.standard_name = "eastward_surface_stokes_drift" ;
        ucur.globwave_name = "eastward_surface_stokes_drift" ;
        ucur.units = "m s-1" ;
        #ucur.scale_factor = 1.f ;
        #ucur.add_offset = 0.f ;
        #ucur.valid_min = -990 ;
        #ucur.valid_max = 990 ;
        ucur.comment = "cur=sqrt(U**2+V**2)" ;

        vcur = self.ncfile.createVariable('vuss', float32, ('time', 'latitude', 'longitude',),fill_value=9.96921e+36)
        vcur.long_name = "northward surface stokes drift" ;
        vcur.standard_name = "northward_surface_stokes_drift" ;
        vcur.globwave_name = "northward_surface_stokes_drift" ;
        vcur.units = "m s-1" ;
        #ucur.scale_factor = 1.f ;
        #ucur.add_offset = 0.f ;
        #ucur.valid_min = -990 ;
        #ucur.valid_max = 990 ;
        vcur.comment = "cur=sqrt(U**2+V**2)" ;

        time_index=0
        for time in self.coverage.read_axis_t():
            cur = self.coverage.read_variable_surface_stokes_drift_at_time(time)
            ucur[time_index:time_index+1,:,:] = resample_2d_to_grid(self.coverage.read_axis_x(),self.coverage.read_axis_y(),self.lon_reg,self.lat_reg,cur[0])
            vcur[time_index:time_index+1,:,:] = resample_2d_to_grid(self.coverage.read_axis_x(),self.coverage.read_axis_y(),self.lon_reg,self.lat_reg,cur[1])
            time_index += 1


    def resample_variable_wetmask(self):
        if self.ncfile == None:
            raise IOError("Please call write_axis() first")

        logging.info('[CoverageInterpolator] Resample variable \'wetmask\' at resolution '+str(self.targetResX)+'/'+str(self.targetResY)+'.')

        var = self.ncfile.createVariable('wetmask', float32, ('time', 'latitude', 'longitude',),fill_value=9.96921e+36)
        var.long_name = "wet_mask" ;
        var.standard_name = "wet_mask" ;
        var.globwave_name = "wet_mask" ;
        var.units = "m" ;
        #wlv.scale_factor = "1.f" ;
        #wlv.add_offset = "0.f" ;
        #wlv.valid_min = "0f" ;
        #wlv.valid_max = 10000f ;

        time_index=0
        for time in self.coverage.read_axis_t():
            var[time_index:time_index+1,:,:] = resample_2d_to_grid(self.coverage.read_axis_x(),self.coverage.read_axis_y(),self.lon_reg,self.lat_reg,self.coverage.read_variable_wetmask_at_time(time))
            time_index += 1


    def resample_variable_sea_water_salinity(self):
        """
        Interpole un champ de courant au niveau donnée à la construction
        @vertical_method: méthode d'interpolation verticale. "nearest" ou "linear"
        """
        if self.ncfile == None:
            raise IOError("Please call write_axis() first")

        var = self.ncfile.createVariable('salinity', float32, ('time', 'depth', 'latitude', 'longitude',),fill_value=9.96921e+36)
        var.long_name = "salinity" ;
        var.standard_name = "salinity" ;
        var.globwave_name = "salinity" ;
        var.units = "psu" ;
        #ucur.scale_factor = 1.f ;
        #ucur.add_offset = 0.f ;
        #ucur.valid_min = -990 ;
        #ucur.valid_max = 990 ;

        time_index=0
        for time in self.coverage.read_axis_t():

            logging.info('[CoverageInterpolator] Resample variable \'salinity\' at time '+str(time)+' to the resolution '+str(self.targetResX)+'/'+str(self.targetResY)+'.')

            level_index = 0
            for level in self.targetDepths:

                if type(level) == int:
                    logging.info('[CoverageInterpolator] At index level '+str(level)+'')
                else:
                    logging.info('[CoverageInterpolator] At depth '+str(level)+' m.')

                var[time_index:time_index+1,level_index:level_index+1,:,:] = resample_2d_to_grid(self.coverage.read_axis_x(),self.coverage.read_axis_y(),self.lon_reg,self.lat_reg,self.coverage.read_variable_salinity_at_time_and_depth(time,level))

                level_index+= 1


            time_index += 1

    def resample_variable_sea_water_temperature(self):
        """
        Interpole un champ de courant au niveau donnée à la construction
        @vertical_method: méthode d'interpolation verticale. "nearest" ou "linear"
        """
        if self.ncfile == None:
            raise IOError("Please call write_axis() first")

        var = self.ncfile.createVariable('temperature', float32, ('time', 'depth', 'latitude', 'longitude',),
                                         fill_value=9.96921e+36)
        var.long_name = "temperature";
        var.standard_name = "temperature";
        var.globwave_name = "temperature";
        var.units = "deg";
        # ucur.scale_factor = 1.f ;
        # ucur.add_offset = 0.f ;
        # ucur.valid_min = -990 ;
        # ucur.valid_max = 990 ;

        time_index = 0
        for time in self.coverage.read_axis_t():

            logging.info('[CoverageInterpolator] Resample variable \'temperature\' at time ' + str(
                time) + ' to the resolution ' + str(self.targetResX) + '/' + str(self.targetResY) + '.')

            level_index = 0
            for level in self.targetDepths:

                if type(level) == int:
                    logging.info('[CoverageInterpolator] At index level ' + str(level) + '')
                else:
                    logging.info('[CoverageInterpolator] At depth ' + str(level) + ' m.')

                var[time_index:time_index + 1, level_index:level_index + 1, :, :] = resample_2d_to_grid(
                    self.coverage.read_axis_x(), self.coverage.read_axis_y(), self.lon_reg, self.lat_reg,
                    self.coverage.read_variable_temperature_at_time_and_depth(time, level))

                level_index += 1

            time_index += 1

    def resample_variable_wind_stress(self):

            if self.ncfile == None:
                raise IOError("Please call write_axis() first")

            logging.info('[CoverageInterpolator] Resample variable \'wind stress\' to the resolution '+str(self.targetResX)+'/'+str(self.targetResY)+'.')

            ucur = self.ncfile.createVariable('uwind_stress', float32, ('time', 'latitude', 'longitude',),fill_value=9.96921e+36)
            ucur.long_name = "eastward wind_stress" ;
            ucur.standard_name = "eastward_wind_stress_velocity" ;
            ucur.globwave_name = "eastward_wind_stress_velocity" ;
            ucur.units = "m s-1" ;
            #ucur.scale_factor = 1.f ;
            #ucur.add_offset = 0.f ;
            #ucur.valid_min = -990 ;
            #ucur.valid_max = 990 ;
            #ucur.comment = "cur=sqrt(U**2+V**2)" ;

            vcur = self.ncfile.createVariable('vwind_stress', float32, ('time', 'latitude', 'longitude',),fill_value=9.96921e+36)
            vcur.long_name = "northward wind_stress" ;
            vcur.standard_name = "northward_wind_stress_velocity" ;
            vcur.globwave_name = "northward_wind_stress_velocity" ;
            vcur.units = "m s-1" ;
            #ucur.scale_factor = 1.f ;
            #ucur.add_offset = 0.f ;
            #ucur.valid_min = -990 ;
            #ucur.valid_max = 990 ;
            vcur.comment = "cur=sqrt(U**2+V**2)" ;

            time_index=0
            for time in self.coverage.read_axis_t():
                cur = self.coverage.read_variable_wind_stress_at_time(time)
                ucur[time_index:time_index+1,:,:] = resample_2d_to_grid(self.coverage.read_axis_x(),self.coverage.read_axis_y(),self.lon_reg,self.lat_reg,cur[0])
                vcur[time_index:time_index+1,:,:] = resample_2d_to_grid(self.coverage.read_axis_x(),self.coverage.read_axis_y(),self.lon_reg,self.lat_reg,cur[1])
                time_index += 1
