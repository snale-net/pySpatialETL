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
from spatialetl.coverage.io import ECMWFReader as CovReader
from netCDF4 import num2date
from scipy.io import loadmat
import numpy as np
from datetime import datetime,timedelta,timezone
import pytz
from spatialetl.utils.logger import logging
from spatialetl.utils.distance import distance_on_unit_sphere

class ECMWFReader(MultiPointReader):

    def __init__(self,myFile,xy,names=None):
        MultiPointReader.__init__(self, myFile);
        self.reader = CovReader(self.filename)
        self.nbPoints = np.shape(xy)[0]
        self.xy_coords = np.zeros([self.nbPoints,2],dtype=np.int32)
        self.xy_values = np.zeros([self.nbPoints, 2])
        self.meta_data = ""

        if names is None:
            nbPoints = np.shape(self.read_axis_x())[0]
            self.names = np.empty([nbPoints], dtype=object)
            for count in range(0, nbPoints):
                self.names[count] = "Point-" + str(count)
        else:
            self.names = names

        for i in range(0, self.nbPoints):
            nearestPoint = self.find_point_index(xy[i][0], xy[i][1])
            logging.info("Nearest point : " + str(nearestPoint[2]) + " / " + str(nearestPoint[3]) + " at " + str(round(nearestPoint[4], 4)) + " km")
            self.meta_data = self.meta_data + "\n# "+str(self.names[i])+" : nearest point in ECMWF file is "+ str(round(nearestPoint[4], 4)) + " km from the target point"
            logging.info("Nearest point (i,j) : " + str(nearestPoint[0]) + " / " + str(nearestPoint[1]))
            self.xy_coords[i] = [nearestPoint[0],nearestPoint[1]]
            self.xy_values[i] = [nearestPoint[2], nearestPoint[3]]

    def close(self):
        self.reader.close()

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
        lon = self.reader.read_axis_x(0, x_size, 0, y_size)
        lat = self.reader.read_axis_y(0, x_size, 0, y_size)

        try:
            mask = self.reader.read_variable_2D_sea_binary_mask(0, x_size, 0, y_size)
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
    def get_t_size(self):
        return self.reader.get_t_size()

    def read_axis_x(self):
        return self.xy_values[:,0]

    def read_axis_y(self):
        return self.xy_values[:,1]

    def read_axis_t(self,timestamp=0):
        return self.reader.read_axis_t(timestamp)

    def read_metadata(self):
        m = {}
        m["data_source"] = "ECMWF file"
        m["meta_data"] = self.meta_data

        return m

    #Scalar
    def read_variable_point_names(self):
        return self.names

    def read_variable_wind_10m_at_time(self,date):

        data = np.zeros([2,self.nbPoints])
        data[:] = np.nan

        all_data = self.reader.read_variable_wind_10m_at_time(date)

        for index_x in range(0, self.nbPoints):
            # comp U
            data[0,index_x] = all_data[0][self.xy_coords[index_x][1],self.xy_coords[index_x][0]]
            # comp V
            data[1,index_x] = all_data[1][self.xy_coords[index_x][1], self.xy_coords[index_x][0]]

        return data

    def read_variable_surface_air_pressure_at_time(self, date):

        data = np.zeros([self.nbPoints])
        data[:] = np.nan

        all_data = self.reader.read_variable_surface_air_pressure_at_time(date)

        for index_x in range(0, self.nbPoints):
            data[index_x] = all_data[self.xy_coords[index_x][1], self.xy_coords[index_x][0]]

        return data

    def read_variable_rainfall_amount_at_time(self, date):

        data = np.zeros([self.nbPoints])
        data[:] = np.nan

        all_data = self.reader.read_variable_rainfall_amount_at_time(date)

        for index_x in range(0, self.nbPoints):
            data[index_x] = all_data[self.xy_coords[index_x][1], self.xy_coords[index_x][0]]

        return data



