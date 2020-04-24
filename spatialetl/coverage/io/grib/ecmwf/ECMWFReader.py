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
from spatialetl.coverage.io.CoverageReader import CoverageReader
from spatialetl.coverage.TimeCoverage import TimeCoverage
from netCDF4 import Dataset,MFDataset, num2date
import numpy as np
from datetime import datetime
from time import strftime
from spatialetl.utils.logger import logging
import os
import cfgrib

class ECMWFReader (CoverageReader):

    def __init__(self, myFile):
        CoverageReader.__init__(self,myFile);

        self.ds = cfgrib.open_file(self.filename)

        lon = self.ds.variables['longitude'].data
        self.new_lon = np.mod(lon + 180.0, 360.0) - 180.0
        lat = self.ds.variables['latitude'].data
        xx, yy = np.meshgrid(self.new_lon, lat)
        self.new_indexes = np.argsort(xx)
        #print(self.ds.attributes['GRIB_edition'])
        #print(sorted(self.ds.dimensions.items()))
        print(sorted(self.ds.variables))

        # if os.path.isfile(self.filename):
        #     self.ncfile = Dataset(self.filename, 'r')
        # elif os.path.isdir(self.filename):
        #     self.ncfile = MFDataset(os.path.join(self.filename, "*.nc"), 'r')
        # elif self.filename.endswith("*"):
        #     self.ncfile = MFDataset(self.filename + ".nc", 'r')
        # else:
        #     raise ValueError("Unable to decode file " + str(self.filename))

    def close(self):
        self.ds.close()

    def is_regular_grid(self):
        return True

    def get_x_size(self):
        return self.ds.dimensions["longitude"]

    def get_y_size(self):
        return self.ds.dimensions["latitude"]

    def get_t_size(self):
        return self.ds.dimensions["time"]

    def read_axis_x(self,xmin,xmax,ymin,ymax):
        return np.array(sorted(self.new_lon)[xmin:xmax])

    def read_axis_y(self,xmin,xmax,ymin,ymax):
        return self.ds.variables['latitude'].data[ymin:ymax]

    def read_axis_t(self,tmin,tmax,timestamp):
        data = self.ds.variables['time'].data[tmin:tmax]
        result = [datetime.utcfromtimestamp(t) \
                    for t in data]

        if timestamp == 1:
            return [(t - TimeCoverage.TIME_DATUM).total_seconds() \
                    for t in result];
        else:
            return result;

    # Variables

    def read_variable_longitude(self,xmin,xmax,ymin,ymax):
        return self.read_axis_x()

    def read_variable_latitude(self,xmin,xmax,ymin,ymax):
        return self.read_axis_y()

    def read_variable_time(self,tmin,tmax):
        return self.read_axis_t(tmin,tmax,timestamp=0)

    def read_variable_2D_land_binary_mask(self,xmin,xmax,ymin,ymax):
        mask = np.ma.filled(self.ds.variables["lsm"].data[0,ymin:ymax,:], fill_value=np.nan)
        return np.take_along_axis(mask, self.new_indexes[ymin:ymax],axis=1)[:,xmin:xmax]

    def read_variable_2D_sea_binary_mask(self,xmin,xmax,ymin,ymax):
        mask = np.ma.filled(self.ds.variables["lsm"].data[0,ymin:ymax,:], fill_value=np.nan)
        mask += 1.0  # inverse le mask
        mask %= 2  # inverse le mask
        return np.take_along_axis(mask, self.new_indexes[ymin:ymax],axis=1)[:,xmin:xmax]

    def read_variable_3D_sea_binary_mask_at_time(self, index_t,xmin,xmax,ymin,ymax):
        mask = np.ma.filled(self.ds.variables["lsm"].data[int(index_t),int(ymin):int(ymax),:], fill_value=np.nan)
        mask += 1.0  # inverse le mask
        mask %= 2  # inverse le mask
        return np.take_along_axis(mask, self.new_indexes[ymin:ymax],axis=1)[:,xmin:xmax]

    def read_variable_3D_land_binary_mask_at_time(self, index_t,xmin,xmax,ymin,ymax):
        mask = np.ma.filled(self.ds.variables["lsm"].data[int(index_t),int(ymin):int(ymax),:], fill_value=np.nan)
        return np.take_along_axis(mask, self.new_indexes[ymin:ymax],axis=1)[:,xmin:xmax]

    #################
    # METEO
    # 2D
    #################

    def read_variable_rainfall_amount_at_time(self, index_t,xmin,xmax,ymin,ymax):
        data = np.ma.filled(self.ds.variables["tp"].data[int(index_t),int(ymin):int(ymax),:], fill_value=np.nan)
        return np.take_along_axis(data, self.new_indexes[ymin:ymax],axis=1)[:,xmin:xmax]

    #################
    # METEO
    # Sea surface
    #################

    def read_variable_surface_air_pressure_at_time(self, index_t,xmin,xmax,ymin,ymax):
        data = np.ma.filled(self.ds.variables["sp"].data[int(index_t),int(ymin):int(ymax),:], fill_value=np.nan)
        data *= 0.01  # Pa to hPa
        return np.take_along_axis(data, self.new_indexes[ymin:ymax],axis=1)[:,xmin:xmax]

    def read_variable_sea_surface_air_pressure_at_time(self,index_t,xmin,xmax,ymin,ymax):
        data = np.ma.filled(self.ds.variables["msl"].data[int(index_t),int(ymin):int(ymax),:], fill_value=np.nan)
        data *= 0.01  # Pa to hPa
        return np.take_along_axis(data, self.new_indexes[ymin:ymax],axis=1)[:,xmin:xmax]

    def read_variable_surface_downward_sensible_heat_flux_at_time(self, index_t,xmin,xmax,ymin,ymax):
        data = np.ma.filled(self.ds.variables["sshf"].data[int(index_t),int(ymin):int(ymax),:], fill_value=np.nan)
        return np.take_along_axis(data, self.new_indexes[ymin:ymax],axis=1)[:,xmin:xmax]

    def read_variable_surface_downward_latent_heat_flux_at_time(self, index_t,xmin,xmax,ymin,ymax):
        data = np.ma.filled(self.ds.variables["slhf"].data[int(index_t),int(ymin):int(ymax),:], fill_value=np.nan)
        return np.take_along_axis(data, self.new_indexes[ymin:ymax],axis=1)[:,xmin:xmax]

    def read_variable_surface_air_temperature_at_time(self, index_t,xmin,xmax,ymin,ymax):
        data = np.ma.filled(self.ds.variables["t2m"].data[int(index_t),int(ymin):int(ymax),:], fill_value=np.nan)
        return np.take_along_axis(data, self.new_indexes[ymin:ymax],axis=1)[:,xmin:xmax]

    def read_variable_dew_point_temperature_at_time(self, index_t,xmin,xmax,ymin,ymax):
        data = np.ma.filled(self.ds.variables["d2m"].data[int(index_t),int(ymin):int(ymax),:], fill_value=np.nan)
        return np.take_along_axis(data, self.new_indexes[ymin:ymax],axis=1)[:,xmin:xmax]

    def read_variable_surface_downwards_solar_radiation_at_time(self, index_t,xmin,xmax,ymin,ymax):
        data = np.ma.filled(self.ds.variables["ssrd"].data[int(index_t),int(ymin):int(ymax),:], fill_value=np.nan)
        return np.take_along_axis(data, self.new_indexes[ymin:ymax],axis=1)[:,xmin:xmax]

    def read_variable_surface_downwards_thermal_radiation_at_time(self, index_t,xmin,xmax,ymin,ymax):
        data = np.ma.filled(self.ds.variables["strd"].data[int(index_t),int(ymin):int(ymax),:], fill_value=np.nan)
        return np.take_along_axis(data, self.new_indexes[ymin:ymax],axis=1)[:,xmin:xmax]

    def read_variable_surface_solar_radiation_at_time(self, index_t,xmin,xmax,ymin,ymax):
        data = np.ma.filled(self.ds.variables["ssr"].data[int(index_t),int(ymin):int(ymax),:], fill_value=np.nan)
        return np.take_along_axis(data, self.new_indexes[ymin:ymax],axis=1)[:,xmin:xmax]

    def read_variable_surface_thermal_radiation_at_time(self, index_t,xmin,xmax,ymin,ymax):
        data = np.ma.filled(self.ds.variables["str"].data[int(index_t), int(ymin):int(ymax),:],fill_value=np.nan)
        return np.take_along_axis(data, self.new_indexes[ymin:ymax],axis=1)[:,xmin:xmax]

    #################
    # METEO
    # At 10 m
    #################

    def read_variable_wind_10m_at_time(self, index_t,xmin,xmax,ymin,ymax):
        u = np.ma.filled(self.ds.variables["u10"].data[int(index_t), int(ymin):int(ymax),:],fill_value=np.nan)
        v = np.ma.filled(self.ds.variables["v10"].data[int(index_t), int(ymin):int(ymax),:],fill_value=np.nan)
        new_u = np.take_along_axis(u, self.new_indexes[ymin:ymax],axis=1)
        new_v = np.take_along_axis(v, self.new_indexes[ymin:ymax],axis=1)
        return [new_u[:,xmin:xmax],new_v[:,xmin:xmax]]
