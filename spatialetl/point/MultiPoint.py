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
import numpy as np
import math
import pandas
from pandas import DatetimeIndex
from datetime import datetime,timedelta
import logging


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

class MultiPoint():
    """"""

    def __init__(self,myReader):
        self.reader = myReader;

        if np.shape(self.reader.read_axis_x())[0] != np.shape(self.reader.read_axis_y())[0] :
            raise ValueError("Longitude axis and latitude axis don't have the same size.")

        self.nb_points = np.shape(self.reader.read_axis_x())[0]

        self.data_source = "Undefined"
        self.name_station = "Undefined"
        self.x_coord = "Undefined"
        self.y_coord = "Undefined"
        self.vertical_datum = "Undefined"
        self.meta_data = "Undefined"

        # try to fill metadata
        self.read_metadata()

        # Read metadata

    def read_metadata(self):
        """
        Lit la metadonnée du fichier si le lecteur contient une fonction read_metadata()
        Returns
        -------

        """

        if "read_metadata" in dir(self.reader) and self.reader.read_metadata() is not None:
            m = self.reader.read_metadata()

            if 'name_station' in m:
                self.name_station = m['name_station']
            if 'data_source' in m:
                self.data_source = m['data_source']
            if 'x_coord' in m:
                self.x_coord = float(m['x_coord'])
            if 'y_coord' in m:
                self.y_coord = float(m['y_coord'])
            if 'vertical_datum' in m:
                self.vertical_datum = m['vertical_datum']
            if 'meta_data' in m:
                self.meta_data = m['meta_data']
        else:
            logging.info("No medatdata available.")
            logging.debug(str(type(self.reader)) + " don't have implemented the function 'read_metadata()'.")

    # Axis
    def get_nb_points(self):
            return self.nb_points

    def read_axis_x(self):
        return self.reader.read_axis_x()

    def read_axis_y(self):
        return self.reader.read_axis_y()

    def find_point_index(self, target_lon, target_lat, method="classic"):
        """Retourne le point le plus proche du point donné en paramètre.
    @param target_lon: Coordonnée longitude du point
    @param target_lat: Coordonnée latitude du point
    @param method : Méthode de calcul. "Classic" = On parcourt toute la grille à la recherche du plus prêt.
    @return: un tableau contenant
     [0] : l'index x du point le plus proche
     [1] : la coordonnée en longitude du point le plus proche
     [2] : la coordonnée en latitude point le plus proche
     [3] : la distance du point le plus proche en kilomètre."""

        min_dist = 10000000
        lon = self.read_axis_x()
        lat = self.read_axis_y()

        if type(target_lon) == int and type(target_lat) == int:

            if (target_lon != target_lat):
                raise ValueError("X index and Y index have to be the sames. Actually X index = " + str(target_lon)+" Actually Y index = " + str(target_lat))

            if target_lon < 0 or target_lon >= self.get_nb_points():
                raise ValueError("XY index have to range between 0 and " + str(
                    self.get_nb_points() - 1) + ". Actually XY index = " + str(target_lon))

            return [target_lon, target_lat, self.read_axis_x()[target_lon], self.read_axis_y()[target_lat], 0.0]

        nearest_point_index = 0
        nearest_lon = np.nan
        nearest_lat = np.nan

        if method == "classic":
            for point in range(0, self.get_nb_points()):

                dist = distance_on_unit_sphere(target_lon, target_lat, lon[point], lat[point])

                if dist < min_dist:
                    min_dist = dist
                    nearest_lon = lon[point]
                    nearest_lat = lat[point]
                    nearest_point_index = point
        else:
            raise RuntimeError("Method " + str(method) + " is not implemented yet.")

        if nearest_lon == np.nan or nearest_lat == np.nan:
            raise RuntimeError("No nearest point found.")

        return [nearest_point_index, nearest_lon, nearest_lat, min_dist]

    # Scalar
    def read_variable_point_names(self):
        return self.reader.read_variable_point_names()

    def read_variable_time(self):
        """
        Read time for all point
        """
        return self.reader.read_variable_time();

    def read_variable_bathymetry(self):
        """
        Read bathymetry for all point
        """
        return self.reader.read_variable_bathymetry();


    def read_variable_bathymetry_at_location(self,x,y):
        """
        Read bathymetry for a specific point
        """
        index_x = self.find_point_index(x,y)
        return self.reader.read_variable_bathymetry_location(index_x[0]);

    def read_variable_sea_surface_temperature(self):
        """
        Read sea_surface_temperature for all point
        """
        return self.reader.read_variable_sea_surface_temperature();

    def read_variable_sea_surface_salinity(self):

        return  self.reader.read_variable_sea_surface_salinity()

    def read_variable_sea_water_pressure_at_sea_water_surface(self):

        return self.reader.read_variable_sea_water_pressure_at_sea_water_surface()

    def read_variable_sea_surface_density(self):

        return self.reader.read_variable_sea_surface_density()

    def read_variable_sea_water_turbidity(self):

        return self.reader.read_variable_sea_water_turbidity()

    def read_variable_sea_water_electrical_conductivity(self):

        return self.reader.read_variable_sea_water_electrical_conductivity()

