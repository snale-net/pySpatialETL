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
from spatialetl.point.io.MultiPointReader import MultiPointReader
from spatialetl.point.TimeMultiPoint import TimeMultiPoint
from spatialetl.coverage.TimeCoverage import TimeCoverage
from spatialetl.coverage.TimeLevelCoverage import TimeLevelCoverage
from spatialetl.coverage.io.netcdf.symphonie.SymphonieBigDataReader import SymphonieBigDataReader as CovReader
from netCDF4 import num2date
from scipy.io import loadmat
import numpy as np
from datetime import datetime,timedelta,timezone
import pytz
from spatialetl.utils.logger import logging
import math

def distance_on_unit_sphere(long1, lat1, long2, lat2):
    """
    Calcule la distance en kilomètre entre les deux point. Les coordonnées
    sont données en longitude, latitude.

    @type  long1: number
    @param long1: Coordonnée X du point 1.
    @type  lat1: number
    @param lat1: Coordonnée Y du point 1.
    @type  long2: number
    @param long2: Coordonnée X du point 1.
    @type  lat2: number
    @param lat2: Coordonnée Y du point 1.
    @return:  la ditance en kilomètre en ces deux point sur la Terre.
    """

    # Convert latitude and longitude to
    # spherical coordinates in radians.
    degrees_to_radians = math.pi / 180.0

    # phi = 90 - latitude
    phi1 = (90.0 - lat1) * degrees_to_radians
    phi2 = (90.0 - lat2) * degrees_to_radians

    # theta = longitude
    theta1 = long1 * degrees_to_radians
    theta2 = long2 * degrees_to_radians

    # Compute spherical distance from spherical coordinates.

    # For two locations in spherical coordinates
    # (1, theta, phi) and (1, theta', phi')
    # cosine( arc length ) =
    #    sin phi sin phi' cos(theta-theta') + cos phi cos phi'
    # distance = rho * arc length

    cos = (math.sin(phi1) * math.sin(phi2) * math.cos(theta1 - theta2) +
           math.cos(phi1) * math.cos(phi2))
    arc = math.acos(cos)

    # Remember to multiply arc by the radius of the earth
    # in your favorite set of units to get length.
    return arc * 6373


class SymphonieBigDataReader(MultiPointReader):

    def __init__(self,myGrid,myFile,xy,names=None):
        MultiPointReader.__init__(self, myFile);
        self.reader = CovReader(myGrid,self.filename)
        self.nbPoints = np.shape(xy)[0]
        self.xy_coords = np.zeros([self.nbPoints,2],dtype=np.int32)
        self.xy_values = np.zeros([self.nbPoints, 2])
        self.meta_data = ""

        if names is None:
            self.names = np.empty([self.nbPoints], dtype=object)
            for count in range(0, self.nbPoints):
                self.names[count] = "Point-" + str(count)
        else:
            self.names = names

        for i in range(0, np.shape(self.xy_coords)[0]):
            nearestPoint = self.find_point_index(xy[i][0], xy[i][1])
            logging.info("Nearest point : " + str(nearestPoint[2]) + " / " + str(nearestPoint[3]) + " at " + str(round(nearestPoint[4], 4)) + " km")
            self.meta_data = self.meta_data + "\n# " + str(self.names[i]) + " : nearest point in SYMPHONIE file is " + str(
                round(nearestPoint[4], 4)) + " km from the target point"
            logging.info("Nearest point (i,j) : " + str(nearestPoint[0]) + " / " + str(nearestPoint[1]))
            self.xy_coords[i] = [nearestPoint[0],nearestPoint[1]]
            self.xy_values[i] = [nearestPoint[2], nearestPoint[3]]

        self.names = names

    def find_point_index(self, target_lon, target_lat, method="classic", only_mask_value=True):
        """Retourne le point le plus proche du point donné en paramètre.
    @param target_lon: Coordonnée longitude du point
    @param target_lat: Coordonnée latitude du point
    @param method : Méthode de calcul. "Classic" = On parcourt toute la grille à la recherche du plus prêt.
    @return: un tableau contenant
     [0] : l'index x du point le plus proche
     [1] : l'index y du point le plus proche
     [2] : la coordonnée en longitude du point le plus proche
     [3] : la coordonnée en latitude point le plus proche
     [4] : la distance du point le plus proche en kilomètre."""
        x_size = self.reader.get_x_size()
        y_size = self.reader.get_y_size()
        lon = self.reader.read_axis_x(0,x_size,0,y_size)
        lat = self.reader.read_axis_y(0,x_size,0,y_size)

        try:
            mask = self.reader.read_variable_2D_sea_binary_mask(0,x_size,0,y_size)
        except NotImplementedError:
            logging.warning("No 2D sea binary mask found")
            # mask = np.ones([self.source_global_y_size, self.source_global_x_size])
            only_mask_value = False

        dist = np.zeros([y_size, x_size])
        dist[:] = 100000

        if method == "classic":
            for x in range(0, x_size):
                for y in range(0, y_size):

                    if only_mask_value:
                        if (mask[y, x] == 1):  # =Terre
                            if self.reader.is_regular_grid():
                                dist[y, x] = distance_on_unit_sphere(target_lon, target_lat, lon[x], lat[y])
                            else:
                                dist[y, x] = distance_on_unit_sphere(target_lon, target_lat, lon[y, x], lat[y, x])
                    else:
                        if self.reader.is_regular_grid():
                            dist[y, x] = distance_on_unit_sphere(target_lon, target_lat, lon[x], lat[y])
                        else:
                            dist[y, x] = distance_on_unit_sphere(target_lon, target_lat, lon[y, x], lat[y, x])

            nearest_y_index, nearest_x_index = np.where(dist == np.min(dist))

            if len(nearest_y_index) == 0 or len(nearest_x_index) == 0:
                raise RuntimeError("No nearest point found.")

            nearest_x_index = nearest_x_index[0]
            nearest_y_index = nearest_y_index[0]
            min_dist = dist[nearest_y_index, nearest_x_index]

            if self.reader.is_regular_grid():
                nearest_lon = lon[nearest_x_index]
                nearest_lat = lat[nearest_y_index]
            else:
                nearest_lon = lon[nearest_y_index, nearest_x_index]
                nearest_lat = lat[nearest_y_index, nearest_x_index]

            return [nearest_x_index, nearest_y_index, nearest_lon, nearest_lat, min_dist]

        elif method == "quick":

            if self.is_regular_grid():

                # Longitude : on cherche l'index le plus proche
                array = np.asarray(lon)
                nearest_x_index = (np.abs(array - target_lon)).argmin()

                # Latitude : on cherche l'index le plus proche
                array = np.asarray(lat)
                nearest_y_index = (np.abs(array - target_lat)).argmin()

                nearest_lon = lon[nearest_x_index]
                nearest_lat = lat[nearest_y_index]

                min_dist = distance_on_unit_sphere(target_lon, target_lat, nearest_lon, nearest_lat)

                return [nearest_x_index, nearest_y_index, nearest_lon, nearest_lat, min_dist]

            else:
                raise NotImplementedError("Method " + str(method) + " is not implemented for regular grid.")

        else:
            raise RuntimeError("Method " + str(method) + " is not implemented yet.")

    # Axis
    def get_z_size(self):
        return self.reader.get_z_size()

    def read_axis_x(self):
        return self.xy_values[:,0]

    def read_axis_y(self):
        return self.xy_values[:,1]

    def read_axis_z(self):
        data = np.zeros([self.nbPoints,np.shape(self.reader.read_axis_z())[0]])
        data[:] = np.nan

        all_data = self.reader.read_axis_z()

        for index_x in range(0, self.nbPoints):
            data[index_x] = all_data[:,self.xy_coords[index_x][1], self.xy_coords[index_x][0]]

        return data

    def read_axis_t(self,timestamp=0):
        return self.reader.read_axis_t(0,self.reader.get_t_size(),timestamp)

    def read_metadata(self):
        m = {}
        m["data_source"] = "SYMPHONIE file"
        m["meta_data"] = self.meta_data

        return m

    def read_variable_longitude(self):
        return self.read_axis_x()

    def read_variable_latitude(self):
        return self.read_axis_y()

    def read_variable_depth(self):
        return self.read_axis_z()

    def read_variable_time(self):
        return self.read_axis_t(timestamp=0)

    def read_variable_point_names(self):
        return self.names

    #################
    # HYDRO
    # Sea Surface
    #################
    def read_variable_sea_surface_height_above_mean_sea_level_at_time(self,index_t):

        data = np.zeros([self.nbPoints])
        data[:] = np.nan

        for index_x in range(0, self.nbPoints):
            data[index_x] = self.reader.read_variable_sea_surface_height_above_mean_sea_level_at_time(index_t,self.xy_coords[index_x][0],self.xy_coords[index_x][0]+1,self.xy_coords[index_x][1],self.xy_coords[index_x][1]+1)

        return data

    def read_variable_sea_surface_temperature_at_time(self,index_t):
        data = np.zeros([self.nbPoints])
        data[:] = np.nan

        for index_x in range(0, self.nbPoints):
            data[index_x] = self.reader.read_variable_sea_surface_temperature_at_time(index_t,self.xy_coords[index_x][0],self.xy_coords[index_x][0]+1,self.xy_coords[index_x][1],self.xy_coords[index_x][1]+1)

        return data

    def read_variable_sea_surface_salinity_at_time(self,index_t):
        data = np.zeros([self.nbPoints])
        data[:] = np.nan

        for index_x in range(0, self.nbPoints):
            data[index_x] = self.reader.read_variable_sea_surface_salinity_at_time(index_t,self.xy_coords[index_x][0],self.xy_coords[index_x][0]+1,self.xy_coords[index_x][1],self.xy_coords[index_x][1]+1)

        return data

    def read_variable_sea_water_velocity_at_sea_water_surface_at_time(self,index_t):

        data = np.zeros([2, self.nbPoints])
        data[:] = np.nan

        for index_x in range(0, self.nbPoints):
            data[0][index_x] = self.reader.read_variable_sea_water_velocity_at_sea_water_surface_at_time(index_t,self.xy_coords[index_x][0],self.xy_coords[index_x][0]+1,self.xy_coords[index_x][1],self.xy_coords[index_x][1]+1)[0]
            data[1][index_x] = self.reader.read_variable_sea_water_velocity_at_sea_water_surface_at_time(index_t,self.xy_coords[index_x][0],self.xy_coords[index_x][0]+1,self.xy_coords[index_x][1],self.xy_coords[index_x][1]+1)[1]

        return data

    #################
    # HYDRO
    # Ground level
    #################
    def read_variable_sea_water_temperature_at_ground_level_at_time(self,index_t):
        data = np.zeros([self.nbPoints])
        data[:] = np.nan

        for index_x in range(0, self.nbPoints):
            data[index_x] = self.reader.read_variable_sea_water_temperature_at_ground_level_at_time(index_t,self.xy_coords[index_x][0],self.xy_coords[index_x][0]+1,self.xy_coords[index_x][1],self.xy_coords[index_x][1]+1)

        return data
    def read_variable_sea_water_salinity_at_ground_level_at_time(self,index_t):
        data = np.zeros([self.nbPoints])
        data[:] = np.nan

        for index_x in range(0, self.nbPoints):
            data[index_x] = self.reader.read_variable_sea_water_salinity_at_ground_level_at_time(index_t,self.xy_coords[index_x][0],self.xy_coords[index_x][0]+1,self.xy_coords[index_x][1],self.xy_coords[index_x][1]+1)

        return data

    def read_variable_sea_water_velocity_at_ground_level_at_time(self,index_t):

        data = np.zeros([2, self.nbPoints])
        data[:] = np.nan

        for index_x in range(0, self.nbPoints):
            data[0][index_x] = self.reader.read_variable_sea_water_velocity_at_ground_level_at_time(index_t,self.xy_coords[index_x][0],self.xy_coords[index_x][0]+1,self.xy_coords[index_x][1],self.xy_coords[index_x][1]+1)[0]
            data[1][index_x] = self.reader.read_variable_sea_water_velocity_at_ground_level_at_time(index_t,self.xy_coords[index_x][0],self.xy_coords[index_x][0]+1,self.xy_coords[index_x][1],self.xy_coords[index_x][1]+1)[1]

        return data

    #################
    # HYDRO
    # 2D
    #################
    def read_variable_bathymetry(self):
        data = np.zeros([self.nbPoints])
        data[:] = np.nan

        for index_x in range(0, self.nbPoints):
            data[index_x] = self.reader.read_variable_bathymetry(self.xy_coords[index_x][0],self.xy_coords[index_x][0]+1,self.xy_coords[index_x][1],self.xy_coords[index_x][1]+1)

        return data

    def read_variable_barotropic_sea_water_velocity_at_time(self,index_t):
        data = np.zeros([2, self.nbPoints])
        data[:] = np.nan

        for index_x in range(0, self.nbPoints):
            data[0][index_x] = self.reader.read_variable_barotropic_sea_water_velocity_at_time(index_t,self.xy_coords[index_x][0],self.xy_coords[index_x][0]+1,self.xy_coords[index_x][1],self.xy_coords[index_x][1]+1)[0]
            data[1][index_x] = self.reader.read_variable_barotropic_sea_water_velocity_at_time(index_t,self.xy_coords[index_x][0],self.xy_coords[index_x][0]+1,self.xy_coords[index_x][1],self.xy_coords[index_x][1]+1)[1]

        return data

    #################
    # HYDRO
    # 3D
    #################
    def read_variable_sea_water_temperature_at_time_and_depth(self,index_t,index_z):
        data = np.zeros([self.nbPoints])
        data[:] = np.nan

        for index_x in range(0, self.nbPoints):
            data[index_x] = self.reader.read_variable_sea_water_temperature_at_time_and_depth(index_t,index_z,self.xy_coords[index_x][0],self.xy_coords[index_x][0]+1,self.xy_coords[index_x][1],self.xy_coords[index_x][1]+1)

        return data

    def read_variable_sea_water_salinity_at_time_and_depth(self,index_t,index_z):
        data = np.zeros([self.nbPoints])
        data[:] = np.nan

        for index_x in range(0, self.nbPoints):
            data[index_x] = self.reader.read_variable_sea_water_salinity_at_time_and_depth(index_t, index_z,self.xy_coords[index_x][0],self.xy_coords[index_x][0]+1,self.xy_coords[index_x][1],self.xy_coords[index_x][1]+1)

        return data

    def read_variable_baroclinic_sea_water_velocity_at_time_and_depth(self,index_t,index_z):
        data = np.zeros([2, self.nbPoints])
        data[:] = np.nan

        for index_x in range(0, self.nbPoints):
            data[0][index_x] = self.reader.read_variable_baroclinic_sea_water_velocity_at_time_and_depth(index_t, index_z,self.xy_coords[index_x][0],self.xy_coords[index_x][0]+1,self.xy_coords[index_x][1],self.xy_coords[index_x][1]+1)[0]
            data[1][index_x] = self.reader.read_variable_baroclinic_sea_water_velocity_at_time_and_depth(index_t, index_z,self.xy_coords[index_x][0],self.xy_coords[index_x][0]+1,self.xy_coords[index_x][1],self.xy_coords[index_x][1]+1)[1]

        return data

    #################
    # WAVES
    # Sea Surface
    #################
    def read_variable_sea_surface_wave_significant_height_at_time(self,index_t):
        data = np.zeros([self.nbPoints])
        data[:] = np.nan

        for index_x in range(0, self.nbPoints):
            data[index_x] = self.reader.read_variable_sea_surface_wave_significant_height_at_time(index_t,self.xy_coords[index_x][0],self.xy_coords[index_x][0]+1,self.xy_coords[index_x][1],self.xy_coords[index_x][1]+1)

        return data

    def read_variable_sea_surface_wave_mean_period_at_time(self,index_t):
        data = np.zeros([self.nbPoints])
        data[:] = np.nan

        for index_x in range(0, self.nbPoints):
            data[index_x] = self.reader.read_variable_sea_surface_wave_mean_period_at_time(index_t,self.xy_coords[index_x][0],self.xy_coords[index_x][0]+1,self.xy_coords[index_x][1],self.xy_coords[index_x][1]+1)

        return data

    def read_variable_sea_surface_wave_to_direction_at_time(self,index_t):
        data = np.zeros([self.nbPoints])
        data[:] = np.nan

        for index_x in range(0, self.nbPoints):
            data[index_x] = self.reader.read_variable_sea_surface_wave_to_direction_at_time(index_t,self.xy_coords[index_x][0],self.xy_coords[index_x][0]+1,self.xy_coords[index_x][1],self.xy_coords[index_x][1]+1)

        return data

    def read_variable_sea_surface_wave_stokes_drift_velocity_at_time(self,index_t):
        data = np.zeros([2, self.nbPoints])
        data[:] = np.nan

        for index_x in range(0, self.nbPoints):
            data[0][index_x] = self.reader.read_variable_sea_surface_wave_stokes_drift_velocity_at_time(index_t,self.xy_coords[index_x][0],self.xy_coords[index_x][0]+1,self.xy_coords[index_x][1],self.xy_coords[index_x][1]+1)[0]
            data[1][index_x] = self.reader.read_variable_sea_surface_wave_stokes_drift_velocity_at_time(index_t,self.xy_coords[index_x][0],self.xy_coords[index_x][0]+1,self.xy_coords[index_x][1],self.xy_coords[index_x][1]+1)[1]

        return data

    #################
    # WAVES
    # Momentum flux
    #################
    def read_variable_atmosphere_momentum_flux_to_waves_at_time(self,index_t):
        data = np.zeros([2, self.nbPoints])
        data[:] = np.nan

        for index_x in range(0, self.nbPoints):
            data[0][index_x] = self.reader.read_variable_atmosphere_momentum_flux_to_waves_at_time(index_t,self.xy_coords[index_x][0],self.xy_coords[index_x][0]+1,self.xy_coords[index_x][1],self.xy_coords[index_x][1]+1)[0]
            data[1][index_x] = self.reader.read_variable_atmosphere_momentum_flux_to_waves_at_time(index_t,self.xy_coords[index_x][0],self.xy_coords[index_x][0]+1,self.xy_coords[index_x][1],self.xy_coords[index_x][1]+1)[1]

        return data

    def read_variable_waves_momentum_flux_to_ocean_at_time(self,index_t):
        data = np.zeros([2, self.nbPoints])
        data[:] = np.nan

        for index_x in range(0, self.nbPoints):
            data[0][index_x] = self.reader.read_variable_waves_momentum_flux_to_ocean_at_time(index_t,self.xy_coords[index_x][0],self.xy_coords[index_x][0]+1,self.xy_coords[index_x][1],self.xy_coords[index_x][1]+1)[0]
            data[1][index_x] = self.reader.read_variable_waves_momentum_flux_to_ocean_at_time(index_t,self.xy_coords[index_x][0],self.xy_coords[index_x][0]+1,self.xy_coords[index_x][1],self.xy_coords[index_x][1]+1)[1]

        return data

    #################
    # METEO
    # Surface air
    #################
    def read_variable_wind_stress_at_time(self,index_t):
        data = np.zeros([2, self.nbPoints])
        data[:] = np.nan

        for index_x in range(0, self.nbPoints):
            data[0][index_x] = self.reader.read_variable_wind_stress_at_time(index_t,self.xy_coords[index_x][0],self.xy_coords[index_x][0]+1,self.xy_coords[index_x][1],self.xy_coords[index_x][1]+1)[0]
            data[1][index_x] = self.reader.read_variable_wind_stress_at_time(index_t,self.xy_coords[index_x][0],self.xy_coords[index_x][0]+1,self.xy_coords[index_x][1],self.xy_coords[index_x][1]+1)[1]

        return data

    #################
    # METEO
    # At 10 m
    #################
    def read_variable_wind_10m_at_time(self,index_t):
        data = np.zeros([2, self.nbPoints])
        data[:] = np.nan

        for index_x in range(0, self.nbPoints):
            data[0][index_x] = self.reader.read_variable_wind_10m_at_time(index_t,self.xy_coords[index_x][0],self.xy_coords[index_x][0]+1,self.xy_coords[index_x][1],self.xy_coords[index_x][1]+1)[0]
            data[1][index_x] = self.reader.read_variable_wind_10m_at_time(index_t,self.xy_coords[index_x][0],self.xy_coords[index_x][0]+1,self.xy_coords[index_x][1],self.xy_coords[index_x][1]+1)[1]

        return data






