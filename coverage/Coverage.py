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
from mpi4py import MPI
from array_split import shape_split
from scipy.spatial.distance import cdist
from coverage.operator.interpolator.InterpolatorCore import resample_2d_to_grid
import logging
 
def distance_on_unit_sphere(long1, lat1,long2, lat2):
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
    degrees_to_radians = math.pi/180.0
         
    # phi = 90 - latitude
    phi1 = (90.0 - lat1)*degrees_to_radians
    phi2 = (90.0 - lat2)*degrees_to_radians
         
    # theta = longitude
    theta1 = long1*degrees_to_radians
    theta2 = long2*degrees_to_radians
         
    # Compute spherical distance from spherical coordinates.
         
    # For two locations in spherical coordinates 
    # (1, theta, phi) and (1, theta', phi')
    # cosine( arc length ) = 
    #    sin phi sin phi' cos(theta-theta') + cos phi cos phi'
    # distance = rho * arc length
     
    cos = (math.sin(phi1)*math.sin(phi2)*math.cos(theta1 - theta2) + 
           math.cos(phi1)*math.cos(phi2))
    arc = math.acos( cos )
 
    # Remember to multiply arc by the radius of the earth 
    # in your favorite set of units to get length.
    return arc*6373

class Coverage(object):
    """
La classe Coverage représente une couverture spatiale sur l'horizontale. Les point qui représentent cette couverture
peuvent être alignés sur une maille régulière (x,y) ou sur une maille non-régulière ((x1,y1),(x2,y2)). En fonction du
type de maille, les fonctions de lecture des axes retourneront des tableaux à une ou deux dimensions. Pour éviter
un chargement en mémoire de la totalité du fichier, la coverage contient un pointeur vers un lecteur. Les couches
sont donc lues à la demande dans le fichier.

Attention, les axes sont toujours inversés dans les tableaux à cause de NetCDF.
Soit l'axe y en premier puis l'axe x. Exemple : [y,x]

@param  myReader: lecteur de fichier
"""

    HORIZONTAL_INTERPOLATION_METHOD = "linear"
    HORIZONTAL_OVERLAPING_SIZE = 5

    def __init__(self, myReader,bbox=None,resolution_x=None,resolution_y=None):
            
        self.reader = myReader;
        self.source_regular_grid = self.reader.is_regular_grid()
        self.target_regular_grid = self.source_regular_grid
        self.map_mpi = None

        self.horizontal_resampling = False
        self.global_source_x_size = self.reader.get_x_size()
        self.global_source_y_size = self.reader.get_y_size()
        self.global_target_res_x = None
        self.global_target_res_y = None
        self.global_target_axis_x = None
        self.global_target_axis_y = None
        self.global_target_x_size = None
        self.global_target_y_size = None

        self.comm = MPI.COMM_WORLD
        self.size = self.comm.Get_size()
        self.rank = self.comm.Get_rank()

        if resolution_x is not None and resolution_y is not None:

            self.horizontal_resampling = True
            self.target_regular_grid = True

            res = np.mean([resolution_x, resolution_y])
            self.global_target_res_x = res
            self.global_target_res_y = res

            if bbox == None:
                # we compute the destination grid
                Ymin = np.min(self.reader.read_axis_y(0, self.global_source_x_size, 0,
                                                               self.get_global_y_size()))
                Ymax = np.max(self.reader.read_axis_y(0, self.global_source_x_size, 0,
                                                              self.global_source_y_size))
                Xmin = np.min(self.reader.read_axis_x(0, self.global_source_x_size, 0,
                                                               self.global_source_y_size))
                Xmax = np.max(self.reader.read_axis_x(0, self.global_source_x_size, 0,
                                                               self.global_source_y_size))
            else:
                Ymin = bbox[2]
                Ymax = bbox[3]
                Xmin = bbox[0]
                Xmax = bbox[1]

            self.global_target_axis_x = np.arange(Xmin, Xmax, res)
            self.global_target_axis_y = np.arange(Ymin, Ymax, res)

            self.global_target_x_size=len(self.global_target_axis_x)
            self.global_target_y_size=len(self.global_target_axis_y)

            if self.rank == 0:
                logging.info('[horizontal_interpolation] Source grid size : (' + str(self.global_source_x_size) + ", " + str(
                    self.global_source_y_size) + ")")
                logging.info('[horizontal_interpolation] Target grid size : (' + str(self.global_target_x_size) + ", " + str(
                    self.global_target_y_size) + ")")

        else:
            self.global_target_x_size = self.global_source_x_size
            self.global_target_y_size = self.global_source_y_size

        self.create_mpi_map()

        if self.rank==0:
            logging.debug("MPI map:")
        for key in self.map_mpi[self.rank]:
            logging.debug("Proc n°"+str(self.rank)+" "+str(key)+"="+str(self.map_mpi[self.rank][key]))
        logging.debug("---------")

        # try to fill metadata
        self.read_metadata()

    def create_mpi_map(self):

        self.map_mpi = np.empty(self.size,dtype=object)
        source_sample = np.zeros([self.global_source_y_size,self.global_source_x_size])
        target_sample = np.zeros([self.global_target_y_size, self.global_target_x_size])

        # Découpage sur 2 axes
        source_slices = shape_split(source_sample.shape, self.size, axis=[0, 0])

        slice_index = 0
        for slyce in source_slices.flatten():
            slice = tuple(slyce)

            m = np.s_[0:2, 3:4]
            print(m)

            map = {}
            # Grille source
            map["source_global_x_min"] = slice[1].start
            map["source_global_x_max"] = slice[1].stop
            map["source_global_y_min"] = slice[0].start
            map["source_global_y_max"] = slice[0].stop
            map["source_x_size"] = map["source_global_x_max"] - map["source_global_x_min"]
            map["source_y_size"] = map["source_global_y_max"] - map["source_global_y_min"]

            map["source_global_x_min_overlap"] = max(0,map["source_global_x_min"]-Coverage.HORIZONTAL_OVERLAPING_SIZE)
            map["source_global_x_max_overlap"] = min(self.global_source_x_size,map["source_global_x_max"]+Coverage.HORIZONTAL_OVERLAPING_SIZE)
            map["source_global_y_min_overlap"] = max(0,map["source_global_y_min"]-Coverage.HORIZONTAL_OVERLAPING_SIZE)
            map["source_global_y_max_overlap"] = min(self.global_source_y_size,map["source_global_y_max"]+Coverage.HORIZONTAL_OVERLAPING_SIZE)
            map["source_x_size_overlap"] = map["source_global_x_max_overlap"] - map["source_global_x_min_overlap"]
            map["source_y_size_overlap"] = map["source_global_y_max_overlap"] -  map["source_global_y_min_overlap"]

            map["source_x_min"] = Coverage.HORIZONTAL_OVERLAPING_SIZE
            map["source_x_max"] = map["source_x_size_overlap"]-Coverage.HORIZONTAL_OVERLAPING_SIZE
            map["source_y_min"] = Coverage.HORIZONTAL_OVERLAPING_SIZE
            map["source_y_max"] = map["source_y_size_overlap"]-Coverage.HORIZONTAL_OVERLAPING_SIZE

            if  map["source_global_x_min"] == 0 :
                map["source_x_min"] = 0

            if  map["source_global_y_min"] == 0 :
                map["source_y_min"] = 0

            if map["source_global_x_max"] == self.global_source_x_size:
                map["source_x_max"] = map["source_x_size_overlap"]

            if map["source_global_y_max"] == self.global_source_y_size:
                 map["source_y_max"] = map["source_y_size_overlap"]

            # Grille destination
            lon = self.reader.read_axis_y(map["source_global_x_min"], map["source_global_x_max"], map["source_global_y_min"],
                                                  map["source_global_y_max"])
            lat = self.reader.read_axis_y(map["source_global_x_min"], map["source_global_x_max"], map["source_global_y_min"],
                                                  map["source_global_y_max"])
            Ymin = np.min(lat)
            Ymax = np.max(lat)
            Xmin = np.min(lon)
            Xmax = np.max(lon)

            lBconer = self.find_point_index(Xmin,Ymin)
            hRcorner = self.find_point_index(Xmax,Ymax)

            map["target_global_x_min"] = lBconer[0]
            map["target_global_x_max"] = hRcorner[0]
            map["target_global_y_min"] = lBconer[1]
            map["target_global_y_max"] = hRcorner[1]

            self.map_mpi[slice_index] = map

            slice_index = slice_index + 1

    # Read metadata
    def read_metadata(self):
        """
        Lit la metadonnée du fichier si le lecteur contient une fonction read_metadata()
        Returns
        -------

        """

        if  "read_metadata" in dir(self.reader):
            m = self.reader.read_metadata()

    def get_x_size(self):
        return self.map_mpi[self.rank]["x_size"]

    def get_y_size(self):
        return self.map_mpi[self.rank]["y_size"]

    def get_global_x_size(self):
        return self.global_target_x_size

    def get_global_y_size(self):
        return self.global_target_y_size

    def is_regular_grid(self):
        """Retourne vrai si la maille est régulière, sinon faux.
    @return:  vrai si la maille est régulière sinon faux."""

        return self.target_regular_grid
    # Axis        
    def read_axis_x(self,type="target",with_overlap=False):
        """Retourne les valeurs (souvent la longitude) de l'axe x.
    @return:  un tableau à une ou deux dimensions selon le type de maille des valeurs de l'axe x (souvent la longitude) : [x] ou [y,x]."""

        if (type == "target_global"):
            return self.global_target_axis_x

        if (type=="source" or self.horizontal_resampling == False) or (type=="target" and self.horizontal_resampling == False):

            data = self.reader.read_axis_x(self.map_mpi[self.rank]["source_global_x_min_overlap"],
                                           self.map_mpi[self.rank]["source_global_x_max_overlap"],
                                           self.map_mpi[self.rank]["source_global_y_min_overlap"],
                                           self.map_mpi[self.rank]["source_global_y_max_overlap"])

            if with_overlap is False:
                if self.source_regular_grid:
                    data = data[self.map_mpi[self.rank]["source_x_min"]:self.map_mpi[self.rank]["source_x_max"]]
                else:
                    data = data[self.map_mpi[self.rank]["source_y_min"]:self.map_mpi[self.rank]["source_y_max"],
                           self.map_mpi[self.rank]["source_x_min"]:self.map_mpi[self.rank]["source_x_max"]]

            return data
        else:
            return self.global_target_axis_x[self.map_mpi[self.rank]["target_global_x_min"]:self.map_mpi[self.rank]["target_global_x_max"]]
        
    def read_axis_y(self,type="target",with_overlap=False):
        """Retourne les valeurs (souvent la latitude) de l'axe y.
    @return:  un tableau à une ou deux dimensions selon le type de maille des valeurs de l'axe y (souvent la latitude) : [x] ou [y,x]."""

        if(type=="target_global"):
            return self.global_target_axis_y

        if (type == "source" or self.horizontal_resampling == False) or (
                type == "target" and self.horizontal_resampling == False):

            data = self.reader.read_axis_y(self.map_mpi[self.rank]["source_global_x_min_overlap"],
                                           self.map_mpi[self.rank]["source_global_x_max_overlap"],
                                           self.map_mpi[self.rank]["source_global_y_min_overlap"],
                                           self.map_mpi[self.rank]["source_global_y_max_overlap"])

            if with_overlap is False:
                if self.source_regular_grid:
                    data = data[self.map_mpi[self.rank]["source_y_min"]:self.map_mpi[self.rank]["source_y_max"]]
                else:
                    data = data[self.map_mpi[self.rank]["source_y_min"]:self.map_mpi[self.rank]["source_y_max"],
                           self.map_mpi[self.rank]["source_x_min"]:self.map_mpi[self.rank]["source_x_max"]]

            return data
        else:
            return self.global_target_axis_y[self.map_mpi[self.rank]["target_global_y_min"]:self.map_mpi[self.rank]["target_global_y_max"]]
        
    def find_point_index(self,target_lon,target_lat,method="classic"):
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
        lon = self.read_axis_x(type="target_global")
        lat = self.read_axis_y(type="target_global")
        #mask = self.read_variable_2D_sea_binary_mask()
        dist = np.zeros([self.global_target_y_size,self.global_target_x_size])
        dist[:] = 100000

        if method=="classic":
            for x in range(0, self.global_target_x_size):
                for y in range(0, self.global_target_y_size):

                    #if(mask[y,x] == 1): #=Terre
                    if self.is_regular_grid():
                        dist[y,x] = distance_on_unit_sphere(target_lon,target_lat,lon[x],lat[y])
                    else:
                        dist[y,x] = distance_on_unit_sphere(target_lon,target_lat,lon[y,x],lat[y,x])

            nearest_y_index,nearest_x_index = np.where(dist == np.min(dist))

            if len(nearest_y_index) == 0 or len(nearest_x_index) == 0:
                raise RuntimeError("No nearest point found.")

            nearest_x_index = nearest_x_index[0]
            nearest_y_index = nearest_y_index[0]
            min_dist = dist[nearest_y_index,nearest_x_index]

            if self.is_regular_grid():
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

                min_dist = distance_on_unit_sphere(target_lon, target_lat,nearest_lon,nearest_lat)

                return [nearest_x_index, nearest_y_index, nearest_lon, nearest_lat, min_dist]

            else:
                raise NotImplementedError("Method " + str(method) + " is not implemented for regular grid.")

        else:
            raise RuntimeError("Method "+str(method)+" is not implemented yet.")
    
    # Variables
    #################
    # HYDRO
    # 2D
    #################
    def read_variable_bathymetry(self):     
        """Retourne la bathymétrie sur toute la couverture
    @return: un tableau en deux dimensions [y,x]."""
        return self.reader.read_variable_bathymetry()
    
    def read_variable_mesh_size(self):     
        """Retourne la taille de la grille sur toute la couverture
    @return: un tableau en deux dimensions [y,x]."""
        return self.reader.read_variable_mesh_size()
    
    def read_variable_2D_sea_binary_mask(self):
        """Retourne le masque terre/mer sur toute la couverture
    @return: un tableau en deux dimensions [y,x].
            0 = Terre
            1 = Mer
    """

        data = self.reader.read_variable_2D_sea_binary_mask(
            self.map_mpi[self.rank]["source_global_x_min_overlap"],
            self.map_mpi[self.rank]["source_global_x_max_overlap"],
            self.map_mpi[self.rank]["source_global_y_min_overlap"],
            self.map_mpi[self.rank]["source_global_y_max_overlap"])

        if self.horizontal_resampling:

            data = resample_2d_to_grid(self.read_axis_x(type="source",with_overlap=True),
                                       self.read_axis_y(type="source",with_overlap=True),
                                       self.read_axis_x(type="target",with_overlap=False),
                                       self.read_axis_y(type="target",with_overlap=False),
                                       data,
                                       Coverage.HORIZONTAL_INTERPOLATION_METHOD)


            return data

        else:
            return data[
               self.map_mpi[self.rank]["source_y_min"]:self.map_mpi[self.rank]["source_y_max"],
               self.map_mpi[self.rank]["source_x_min"]:self.map_mpi[self.rank]["source_x_max"]]


    #################
    # METEO
    # 2D
    #################
    def read_variable_topography(self):     
        """Retourne la topographie sur toute la couverture
    @return: un tableau en deux dimensions [y,x]."""
        return self.reader.read_variable_topography()

    #OTHERS
    def read_variable_Ha(self):
        """Retourne l'amplitude de la réanalyse
    @return: un tableau en deux dimensions [y,x]."""
        return self.reader.read_variable_Ha()

    # def read_variable_sea_surface_density(self):
    #     """Retourne la densité de l'eau de surface
    #         @return: un tableau en deux dimensions [y,x]."""
    #     return self.reader.read_variable_sea_surface_density()
    #
    # def read_variable_sea_water_turbidity(self):
    #     """Retourne la turbidité de l'eau de surface
    #                 @return: un tableau en deux dimensions [y,x]."""
    #     return self.reader.read_variable_sea_water_turbidity()

        
        