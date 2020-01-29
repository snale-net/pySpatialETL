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
from netCDF4 import Dataset
from netCDF4 import date2num
from numpy import float32
from numpy import float64
import numpy as np
import logging
from spatialetl.coverage.io import WW3UnstructuredReader

class OASISWriter(object):

     def __init__(self,symphonie,ww3):
        self.symphonieCoverage = symphonie
        self.ww3Coverage = ww3

     def write_grids(self,filename):

        logging.info('[OASISWriter] Write grids.')

        ncfile = Dataset(filename, 'w', format='NETCDF4')
        ncfile.description = 'OASIS Writer. Generated with Coverage Processing tools'

        # dimensions SYMPHONIE
        ncfile.createDimension('y_symt', self.symphonieCoverage.get_y_size())
        ncfile.createDimension('x_symt', self.symphonieCoverage.get_x_size())


        symphonieLatitudes = ncfile.createVariable('symt.lat', float64, ('y_symt','x_symt',))
        symphonieLatitudes.units = "degree_north" ;
        symphonieLatitudes.long_name = "latitude" ;
        symphonieLatitudes.standard_name = "latitude" ;
        symphonieLatitudes.valid_min = "-90.f";
        symphonieLatitudes.valid_max = "90.f" ;
        symphonieLatitudes.axis = "Y" ;

        symphonieLongitudes = ncfile.createVariable('symt.lon', float64, ('y_symt','x_symt',))
        symphonieLongitudes.units = "degree_east" ;
        symphonieLongitudes.long_name = "longitude" ;
        symphonieLongitudes.standard_name = "longitude" ;
        symphonieLongitudes.valid_min = "-180.f" ;
        symphonieLongitudes.valid_max = "180.f" ;
        symphonieLongitudes.axis = "X" ;

        # data SYMPHONIE
        symphonieLatitudes[:,:] = self.symphonieCoverage.read_axis_y();
        symphonieLongitudes[:,:] = self.symphonieCoverage.read_axis_x();

        if self.ww3Coverage.is_regular_grid()==True:

            # dimensions WW3
            ncfile.createDimension('y_ww3t', self.ww3Coverage.get_y_size())
            ncfile.createDimension('x_ww3t', self.ww3Coverage.get_x_size())

            # variables
            ww3Latitudes = ncfile.createVariable('ww3t.lat', float64, ('y_ww3t',))
            ww3Latitudes.units = "degree_north" ;
            ww3Latitudes.long_name = "latitude" ;
            ww3Latitudes.standard_name = "latitude" ;
            ww3Latitudes.valid_min = "-90.f";
            ww3Latitudes.valid_max = "90.f" ;
            ww3Latitudes.axis = "Y" ;

            ww3Longitudes = ncfile.createVariable('ww3t.lon', float64, ('x_ww3t',))
            ww3Longitudes.units = "degree_east" ;
            ww3Longitudes.long_name = "longitude" ;
            ww3Longitudes.standard_name = "longitude" ;
            ww3Longitudes.valid_min = "-180.f" ;
            ww3Longitudes.valid_max = "180.f" ;
            ww3Longitudes.axis = "X" ;

            # data
            ww3Latitudes[:] = self.ww3Coverage.read_axis_y();
            ww3Longitudes[:] = self.ww3Coverage.read_axis_x();

        else:

            if(isinstance(self.ww3Coverage.reader, WW3UnstructuredReader)):

                # dimensions WW3
                ncfile.createDimension('y_ww3t', 1)
                ncfile.createDimension('x_ww3t', self.ww3Coverage.get_x_size()*self.ww3Coverage.get_y_size())

                ww3Latitudes = ncfile.createVariable('ww3t.lat', float64, ('y_ww3t','x_ww3t',))
                ww3Latitudes.units = "degree_north" ;
                ww3Latitudes.long_name = "latitude" ;
                ww3Latitudes.standard_name = "latitude" ;
                ww3Latitudes.valid_min = "-90.f";
                ww3Latitudes.valid_max = "90.f" ;
                ww3Latitudes.axis = "Y" ;

                ww3Longitudes = ncfile.createVariable('ww3t.lon', float64, ('y_ww3t','x_ww3t',))
                ww3Longitudes.units = "degree_east" ;
                ww3Longitudes.long_name = "longitude" ;
                ww3Longitudes.standard_name = "longitude" ;
                ww3Longitudes.valid_min = "-180.f" ;
                ww3Longitudes.valid_max = "180.f" ;
                ww3Longitudes.axis = "X" ;

                # data
                ww3Latitudes[:,:] = np.reshape(self.ww3Coverage.read_axis_y(), (1, self.ww3Coverage.get_x_size()*self.ww3Coverage.get_y_size()));
                ww3Longitudes[:,:] = np.reshape(self.ww3Coverage.read_axis_x(), (1, self.ww3Coverage.get_x_size()*self.ww3Coverage.get_y_size()))

            else:
                # dimensions WW3
                ncfile.createDimension('y_ww3t', self.ww3Coverage.get_y_size())
                ncfile.createDimension('x_ww3t', self.ww3Coverage.get_x_size())

                ww3Latitudes = ncfile.createVariable('ww3t.lat', float64, ('y_ww3t','x_ww3t',))
                ww3Latitudes.units = "degree_north" ;
                ww3Latitudes.long_name = "latitude" ;
                ww3Latitudes.standard_name = "latitude" ;
                ww3Latitudes.valid_min = "-90.f";
                ww3Latitudes.valid_max = "90.f" ;
                ww3Latitudes.axis = "Y" ;

                ww3Longitudes = ncfile.createVariable('ww3t.lon', float64, ('y_ww3t','x_ww3t',))
                ww3Longitudes.units = "degree_east" ;
                ww3Longitudes.long_name = "longitude" ;
                ww3Longitudes.standard_name = "longitude" ;
                ww3Longitudes.valid_min = "-180.f" ;
                ww3Longitudes.valid_max = "180.f" ;
                ww3Longitudes.axis = "X" ;

                # data
                ww3Latitudes[:,:] = self.ww3Coverage.read_axis_y();
                ww3Longitudes[:,:] = self.ww3Coverage.read_axis_x();

        ncfile.close()

     def write_masks(self,filename):

        logging.info('[OASISWriter] Write masks.')

        ncfile = Dataset(filename, 'w', format='NETCDF4')
        ncfile.description = 'OASIS Writer. Generated with Coverage Processing tools'

        # dimensions SYMPHONIE
        ncfile.createDimension('y_symt', self.symphonieCoverage.get_y_size())
        ncfile.createDimension('x_symt', self.symphonieCoverage.get_x_size())


        symphonieMask = ncfile.createVariable('symt.msk', np.int32, ('y_symt','x_symt',))
        symphonieMask.units = "" ;
        symphonieMask.long_name = "Land/Sea Mask" ;
        symphonieMask.standard_name = "land_sea_mask" ;
        symphonieMask.valid_min = "0";
        symphonieMask.valid_max = "1" ;


        # data SYMPHONIE
        data = self.symphonieCoverage.read_variable_2D_mask();
        for i in range(0, self.symphonieCoverage.get_x_size()):
            for j in range(0, self.symphonieCoverage.get_y_size()):
                if (data[j,i]== 0):
                    symphonieMask[j,i] = 1
                else:
                    symphonieMask[j,i] = 0

        if self.ww3Coverage.is_regular_grid()==True:

            # dimensions WW3
            ncfile.createDimension('y_ww3t', self.ww3Coverage.get_y_size())
            ncfile.createDimension('x_ww3t', self.ww3Coverage.get_x_size())

            # variables
            ww3Mask = ncfile.createVariable('ww3t.msk', np.int32, ('y_ww3t','x_ww3t',))
            ww3Mask.units = "" ;
            ww3Mask.long_name = "Land/Sea Mask" ;
            ww3Mask.standard_name = "land_sea_mask" ;
            ww3Mask.valid_min = "0";
            ww3Mask.valid_max = "1" ;

            # data
            data = self.ww3Coverage.read_variable_2D_mask();
            for i in range(0, self.ww3Coverage.get_x_size()):
                for j in range(0, self.ww3Coverage.get_y_size()):
                    if (data[j,i]== 0):
                        ww3Mask[j,i] = 1
                    else:
                        ww3Mask[j,i] = 0

        else:

            if(isinstance(self.ww3Coverage.reader, WW3UnstructuredReader)):

                # dimensions WW3
                ncfile.createDimension('y_ww3t', 1)
                ncfile.createDimension('x_ww3t', self.ww3Coverage.get_x_size()*self.ww3Coverage.get_y_size())

                # variables
                ww3Mask = ncfile.createVariable('ww3t.msk', np.int32, ('y_ww3t','x_ww3t',))
                ww3Mask.units = "" ;
                ww3Mask.long_name = "Land/Sea Mask" ;
                ww3Mask.standard_name = "land_sea_mask" ;
                ww3Mask.valid_min = "0";
                ww3Mask.valid_max = "1" ;

                print("inverse le mask !!!!")

                # data
                ww3Mask[:,:] = np.reshape(self.ww3Coverage.read_variable_2D_mask(), (1, self.ww3Coverage.get_x_size()*self.ww3Coverage.get_y_size()))

            else:
                # dimensions WW3
                ncfile.createDimension('y_ww3t', self.ww3Coverage.get_y_size())
                ncfile.createDimension('x_ww3t', self.ww3Coverage.get_x_size())

                 # variables
                ww3Mask = ncfile.createVariable('ww3t.msk', np.int32, ('y_ww3t','x_ww3t',))
                ww3Mask.units = "" ;
                ww3Mask.long_name = "Land/Sea Mask" ;
                ww3Mask.standard_name = "land_sea_mask" ;
                ww3Mask.valid_min = "0";
                ww3Mask.valid_max = "1" ;

                # data
                data = self.ww3Coverage.read_variable_2D_mask();
                for i in range(0, self.ww3Coverage.get_x_size()):
                    for j in range(0, self.ww3Coverage.get_y_size()):
                        #ww3Mask[j,i] = 0
                        if (data[j,i]== 0):
                            ww3Mask[j,i] = 1
                        else:
                            ww3Mask[j,i] = 0

        ncfile.close()
       
    