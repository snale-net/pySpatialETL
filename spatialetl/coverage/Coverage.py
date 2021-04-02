#! /usr/bin/env python2.7
# -*- coding: utf-8 -*-
#
# pySpatialETL is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# any later version.
#
# pySpatialETL is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# Author : Fabien Rétif - fabien.retif@zoho.com
#
from __future__ import division, print_function, absolute_import

import numpy as np
from array_split import shape_split
from mpi4py import MPI

from spatialetl.exception.NotFoundInRankError import NotFoundInRankError
from spatialetl.operator.interpolator.InterpolatorCore import resample_2d_to_grid
from spatialetl.utils.distance import distance_on_unit_sphere
from spatialetl.utils.logger import logging


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
    HORIZONTAL_OVERLAPING_SIZE = 2

    def __init__(self, myReader,bbox=None,resolution_x=None,resolution_y=None):
        self.reader = myReader;
        # MPI
        self.map_mpi = None
        self.comm = MPI.COMM_WORLD
        self.size = self.comm.Get_size()
        self.rank = self.comm.Get_rank()

        self.source_regular_grid = self.reader.is_regular_grid()
        self.target_regular_grid = self.source_regular_grid
        self.horizontal_resampling = False

        self.source_global_x_size = self.reader.get_x_size()
        self.source_global_y_size = self.reader.get_y_size()
        self.source_global_axis_x = self.reader.read_axis_x(0, self.source_global_x_size, 0,
                                                            self.source_global_y_size)
        self.source_global_axis_y = self.reader.read_axis_y(0, self.source_global_x_size, 0,
                                                            self.source_global_y_size)

        # On réduit en fonction de la bbox
        if bbox is None:
            # we compute the destination grid
            Ymin = np.min(self.source_global_axis_y)
            Ymax = np.max(self.source_global_axis_y)
            Xmin = np.min(self.source_global_axis_x)
            Xmax = np.max(self.source_global_axis_x)
        else:
            if self.check_bbox_validity(bbox) is False:
                raise ValueError("Your Bbox is not valid or is out side the coverage")
            Ymin = bbox[2]
            Ymax = bbox[3]
            Xmin = bbox[0]
            Xmax = bbox[1]

        if self.is_regular_grid(type="source"):

            idx = np.where((self.source_global_axis_x >= Xmin) &
                           (self.source_global_axis_x <= Xmax))

            xmin = np.min(idx[0])
            xmax = np.max(idx[0]) + 1

            idx = np.where((self.source_global_axis_y >= Ymin) &
                           (self.source_global_axis_y <= Ymax))

            ymin = np.min(idx[0])
            ymax = np.max(idx[0]) + 1

            self.target_global_axis_x = self.source_global_axis_x[xmin:xmax]
            self.target_global_x_size = xmax - xmin
            self.target_global_axis_y = self.source_global_axis_y[ymin:ymax]
            self.target_global_y_size = ymax - ymin

        else:

            idx = np.where((self.source_global_axis_x >= Xmin) &
                           (self.source_global_axis_x <= Xmax) &
                           (self.source_global_axis_y >= Ymin) &
                           (self.source_global_axis_y <= Ymax))

            if len(idx[0]) == 0:
                raise ValueError("No values found")

            ymin = np.min(idx[0])
            ymax = np.max(idx[0])+1
            xmin = np.min(idx[1])
            xmax = np.max(idx[1])+1

            self.target_global_axis_x = self.source_global_axis_x[ymin:ymax, xmin:xmax]
            self.target_global_x_size = xmax - xmin
            self.target_global_axis_y = self.source_global_axis_y[ymin:ymax, xmin:xmax]
            self.target_global_y_size = ymax - ymin

        # On calcule la grille de destination
        self.target_global_res_x = None
        self.target_global_res_y = None

        if resolution_x is not None and resolution_y is not None:

            self.horizontal_resampling = True
            self.target_regular_grid = True

            res = np.mean([resolution_x, resolution_y])
            self.target_global_res_x = res
            self.target_global_res_y = res

            self.target_global_axis_x = np.arange(Xmin, Xmax, res)
            self.target_global_axis_y = np.arange(Ymin, Ymax, res)

            self.target_global_x_size=len(self.target_global_axis_x)
            self.target_global_y_size=len(self.target_global_axis_y)

            if type(self) == Coverage:
                if self.rank == 0:
                    logging.info('[horizontal_interpolation] Source grid size : (' + str(self.source_global_x_size) + ", " + str(
                        self.source_global_y_size) + ")")
                    logging.info('[horizontal_interpolation] Target grid size : (' + str(self.target_global_x_size) + ", " + str(
                        self.target_global_y_size) + ")")

        if type(self)== Coverage:
            self.create_mpi_map()
            self.update_mpi_map()
            if self.rank==0:
                logging.debug("MPI map:")
            for key in self.map_mpi[self.rank]:
                logging.debug("Proc n°"+str(self.rank)+" "+str(key)+"="+str(self.map_mpi[self.rank][key]))
            logging.debug("---------")

        # try to fill metadata
        self.read_metadata()

    def check_bbox_validity(self,candidate):
        Ymin = candidate[2]
        Ymax = candidate[3]
        Xmin = candidate[0]
        Xmax = candidate[1]

        if Xmax <= Xmin:
            return False
        if Ymax <= Ymin:
            return False

        if Ymin < np.min(self.source_global_axis_y):
            return False
        if Ymax > np.max(self.source_global_axis_y):
            return False
        if Xmin < np.min(self.source_global_axis_x):
            return False
        if Xmax > np.max(self.source_global_axis_x):
            return False

        return True

    def check_point_is_inside(self, target_lon, target_lat, lon, lat,tolerance=5):

        if np.round(target_lat,decimals=tolerance) < np.min(np.round(lat,decimals=tolerance)):
            return False
        if np.round(target_lat,decimals=tolerance) > np.max(np.round(lat,decimals=tolerance)):
            return False
        if np.round(target_lon,decimals=tolerance) < np.min(np.round(lon,decimals=tolerance)):
            return False
        if np.round(target_lon,decimals=tolerance) > np.max(np.round(lon,decimals=tolerance)):
            return False

        return True

    def create_mpi_map(self):
        self.map_mpi = np.empty(self.size, dtype=object)
        target_sample = (self.target_global_y_size, self.target_global_x_size)

        # Découpage des axes
        target_slices = shape_split(target_sample, self.size, axis=[0, 0])

        slice_index = 0
        for slyce in target_slices.flatten():
            slice = tuple(slyce)

            map = {}
            # Grille source
            map["dst_global_x"] = slice[1]
            map["dst_global_y"] = slice[0]

            map["dst_local_x_size"] = map["dst_global_x"].stop - map["dst_global_x"].start
            map["dst_local_y_size"] = map["dst_global_y"].stop - map["dst_global_y"].start

            dst_global_x_min_overlap = max(0, map["dst_global_x"].start - Coverage.HORIZONTAL_OVERLAPING_SIZE)
            dst_global_x_max_overlap = min(self.target_global_x_size,
                                           map["dst_global_x"].stop + Coverage.HORIZONTAL_OVERLAPING_SIZE)
            map["dst_global_x_overlap"] = np.s_[dst_global_x_min_overlap:dst_global_x_max_overlap]

            dst_global_y_min_overlap = max(0, map["dst_global_y"].start - Coverage.HORIZONTAL_OVERLAPING_SIZE)
            dst_global_y_max_overlap = min(self.target_global_y_size,
                                           map["dst_global_y"].stop + Coverage.HORIZONTAL_OVERLAPING_SIZE)
            map["dst_global_y_overlap"] = np.s_[dst_global_y_min_overlap:dst_global_y_max_overlap]

            map["dst_global_x_size_overlap"] = map["dst_global_x_overlap"].stop - map["dst_global_x_overlap"].start
            map["dst_global_y_size_overlap"] = map["dst_global_y_overlap"].stop - map["dst_global_y_overlap"].start

            dst_x_min = Coverage.HORIZONTAL_OVERLAPING_SIZE
            dst_x_max = map["dst_global_x_size_overlap"] - Coverage.HORIZONTAL_OVERLAPING_SIZE
            dst_y_min = Coverage.HORIZONTAL_OVERLAPING_SIZE
            dst_y_max = map["dst_global_y_size_overlap"] - Coverage.HORIZONTAL_OVERLAPING_SIZE

            if map["dst_global_x"].start == 0:
                dst_x_min = 0

            if map["dst_global_x"].stop == self.target_global_x_size:
                dst_x_max = map["dst_global_x_size_overlap"]

            if map["dst_global_y"].start == 0:
                dst_y_min = 0

            if map["dst_global_y"].stop == self.target_global_y_size:
                dst_y_max = map["dst_global_y_size_overlap"]

            map["dst_local_x"] = np.s_[dst_x_min:dst_x_max]
            map["dst_local_y"] = np.s_[dst_y_min:dst_y_max]

            # Source grille
            map["src_global_x"] = map["dst_global_x"]
            map["src_global_y"] = map["dst_global_y"]

            map["src_global_x_overlap"] = map["dst_global_x_overlap"]
            map["src_global_y_overlap"] = map["dst_global_y_overlap"]

            map["src_local_x"] = map["dst_local_x"]
            map["src_local_y"] = map["dst_local_y"]

            map["src_local_x_size"] = map["dst_local_x_size"]
            map["src_local_y_size"] = map["dst_local_y_size"]

            map["src_local_x_size_overlap"] = map["dst_global_x_size_overlap"]
            map["src_local_y_size_overlap"] = map["dst_global_y_size_overlap"]

            self.map_mpi[slice_index] = map

            slice_index = slice_index + 1

    def update_mpi_map(self):

        if self.is_regular_grid(type="source"):

            idx = np.where((self.source_global_axis_x >= np.min(self.read_axis_x(type="target", with_overlap=False))) &
                           (self.source_global_axis_x <= np.max(self.read_axis_x(type="target", with_overlap=False))))

            xmin = np.min(idx[0])
            xmax = np.max(idx[0]) + 1

            idx = np.where((self.source_global_axis_y >= np.min(self.read_axis_y(type="target", with_overlap=False))) &
                           (self.source_global_axis_y <= np.max(self.read_axis_y(type="target", with_overlap=False))))

            ymin = np.min(idx[0])
            ymax = np.max(idx[0]) + 1

        else:

            idx = np.where(
                (self.source_global_axis_x >= np.min(self.read_axis_x(type="target", with_overlap=False))) &
                (self.source_global_axis_x <= np.max(self.read_axis_x(type="target", with_overlap=False))) &
                (self.source_global_axis_y >= np.min(self.read_axis_y(type="target", with_overlap=False))) &
                (self.source_global_axis_y <= np.max(self.read_axis_y(type="target", with_overlap=False))))

            ymin = np.min(idx[0])
            ymax = np.max(idx[0]) + 1
            xmin = np.min(idx[1])
            xmax = np.max(idx[1]) + 1

        # Version 2
        # SRC GLOBAL
        self.map_mpi[self.rank]["src_global_x"] = np.s_[xmin:xmax]
        self.map_mpi[self.rank]["src_global_x_size"] = xmax - xmin
        self.map_mpi[self.rank]["src_global_y"] = np.s_[ymin:ymax]
        self.map_mpi[self.rank]["src_global_y_size"] = ymax - ymin

        dst_global_x_min_overlap = max(0, self.map_mpi[self.rank][
            "src_global_x"].start - Coverage.HORIZONTAL_OVERLAPING_SIZE)
        dst_global_x_max_overlap = min(self.source_global_x_size,
                                       self.map_mpi[self.rank][
                                           "src_global_x"].stop + Coverage.HORIZONTAL_OVERLAPING_SIZE)
        self.map_mpi[self.rank]["src_global_x_overlap"] = np.s_[
                                                          dst_global_x_min_overlap:dst_global_x_max_overlap]

        dst_global_y_min_overlap = max(0, self.map_mpi[self.rank][
            "src_global_y"].start - Coverage.HORIZONTAL_OVERLAPING_SIZE)
        dst_global_y_max_overlap = min(self.source_global_y_size,
                                       self.map_mpi[self.rank][
                                           "src_global_y"].stop + Coverage.HORIZONTAL_OVERLAPING_SIZE)
        self.map_mpi[self.rank]["src_global_y_overlap"] = np.s_[
                                                          dst_global_y_min_overlap:dst_global_y_max_overlap]

        self.map_mpi[self.rank]["src_global_x_size_overlap"] = self.map_mpi[self.rank][
                                                                   "src_global_x_overlap"].stop - \
                                                               self.map_mpi[self.rank][
                                                                   "src_global_x_overlap"].start
        self.map_mpi[self.rank]["src_global_y_size_overlap"] = self.map_mpi[self.rank][
                                                                   "src_global_y_overlap"].stop - \
                                                               self.map_mpi[self.rank][
                                                                   "src_global_y_overlap"].start

        self.map_mpi[self.rank]["src_local_x_size"] = xmax - xmin
        self.map_mpi[self.rank]["src_local_y_size"] = ymax - ymin
        self.map_mpi[self.rank]["src_local_x"] = np.s_[0:self.map_mpi[self.rank]["src_local_x_size"]]
        self.map_mpi[self.rank]["src_local_y"] = np.s_[0:self.map_mpi[self.rank]["src_local_y_size"]]

        # OVERLAP
        self.map_mpi[self.rank]["src_local_x_size_overlap"] = self.map_mpi[self.rank][
            "src_global_x_size_overlap"]
        self.map_mpi[self.rank]["src_local_y_size_overlap"] = self.map_mpi[self.rank][
            "src_global_y_size_overlap"]

        self.map_mpi[self.rank]["src_local_x_overlap"] = np.s_[
                                                         0:self.map_mpi[self.rank]["src_local_x_size_overlap"]]
        self.map_mpi[self.rank]["src_local_y_overlap"] = np.s_[
                                                         0:self.map_mpi[self.rank]["src_local_y_size_overlap"]]

    # Read metadata
    def read_metadata(self):
        """
        Lit la metadonnée du fichier si le lecteur contient une fonction read_metadata()
        Returns
        -------

        """

        if  "read_metadata" in dir(self.reader):
            m = self.reader.read_metadata()

    def get_x_size(self,type="target",with_overlap=False):
        if type == "target_global":
            return self.target_global_x_size
        elif type == "source_global":
            return self.source_global_x_size
        elif type == "source" and with_overlap is True:
            return self.map_mpi[self.rank]["src_local_x_size_overlap"]
        elif type == "source" and with_overlap is False:
            return self.map_mpi[self.rank]["src_local_x_size"]
        elif type == "target" and with_overlap is True:
            return self.map_mpi[self.rank]["dst_local_x_size_overlap"]
        else:
            return self.map_mpi[self.rank]["dst_local_x_size"]

    def get_y_size(self,type="target",with_overlap=False):
        if type == "target_global":
            return self.target_global_y_size
        elif type == "source_global":
            return self.source_global_y_size
        elif type == "source" and with_overlap is True:
            return self.map_mpi[self.rank]["src_local_y_size_overlap"]
        elif type == "source" and with_overlap is False:
            return self.map_mpi[self.rank]["src_local_y_size"]
        elif type == "target" and with_overlap is True:
            return self.map_mpi[self.rank]["dst_local_y_size_overlap"]
        else:
            return self.map_mpi[self.rank]["dst_local_y_size"]

    def is_regular_grid(self,type="target"):
        """Retourne vrai si la maille est régulière, sinon faux.
    @return:  vrai si la maille est régulière sinon faux."""

        if type == "target":
            return self.target_regular_grid
        else:
            return self.source_regular_grid

    # Axis        
    def read_axis_x(self,type="target",with_overlap=False):
        """Retourne les valeurs (souvent la longitude) de l'axe x.
    @return:  un tableau à une ou deux dimensions selon le type de maille des valeurs de l'axe x (souvent la longitude) : [x] ou [y,x]."""

        if type == "target_global":
            return self.target_global_axis_x

        elif type == "source_global":
            return self.source_global_axis_x

        elif type == "source" and with_overlap is True:
            return self.reader.read_axis_x(self.map_mpi[self.rank]["src_global_x_overlap"].start,
                                           self.map_mpi[self.rank]["src_global_x_overlap"].stop,
                                           self.map_mpi[self.rank]["src_global_y_overlap"].start,
                                           self.map_mpi[self.rank]["src_global_y_overlap"].stop)
        elif type == "source" and with_overlap is False:
            return self.reader.read_axis_x(self.map_mpi[self.rank]["src_global_x"].start,
                                    self.map_mpi[self.rank]["src_global_x"].stop,
                                    self.map_mpi[self.rank]["src_global_y"].start,
                                    self.map_mpi[self.rank]["src_global_y"].stop)

        elif type == "target" and with_overlap is True:

            if self.is_regular_grid():
                return self.target_global_axis_x[self.map_mpi[self.rank]["dst_global_x_overlap"]]
            else:
                return self.target_global_axis_x[self.map_mpi[self.rank]["dst_global_y_overlap"],
                                                 self.map_mpi[self.rank]["dst_global_x_overlap"]]
        else :

            if self.is_regular_grid():
                return self.target_global_axis_x[self.map_mpi[self.rank]["dst_global_x"]]
            else:
                return self.target_global_axis_x[self.map_mpi[self.rank]["dst_global_y"],self.map_mpi[self.rank]["dst_global_x"]]
        
    def read_axis_y(self,type="target",with_overlap=False):
        """Retourne les valeurs (souvent la latitude) de l'axe y.
    @return:  un tableau à une ou deux dimensions selon le type de maille des valeurs de l'axe y (souvent la latitude) : [x] ou [y,x]."""

        if type=="target_global":
            return self.target_global_axis_y

        elif type == "source_global":
            return self.source_global_axis_y

        elif type == "source" and with_overlap is True:
            return self.reader.read_axis_y(self.map_mpi[self.rank]["src_global_x_overlap"].start,
                                           self.map_mpi[self.rank]["src_global_x_overlap"].stop,
                                           self.map_mpi[self.rank]["src_global_y_overlap"].start,
                                           self.map_mpi[self.rank]["src_global_y_overlap"].stop)
        elif type == "source" and with_overlap is False:
            return self.reader.read_axis_y(self.map_mpi[self.rank]["src_global_x"].start,
                                           self.map_mpi[self.rank]["src_global_x"].stop,
                                           self.map_mpi[self.rank]["src_global_y"].start,
                                           self.map_mpi[self.rank]["src_global_y"].stop)

        elif type == "target" and with_overlap is True:

            if self.is_regular_grid():
                return self.target_global_axis_y[self.map_mpi[self.rank]["dst_global_y_overlap"]]
            else:
                return self.target_global_axis_y[self.map_mpi[self.rank]["dst_global_y_overlap"],
                                                 self.map_mpi[self.rank]["dst_global_x_overlap"]]
        else:

            if self.is_regular_grid():
                return self.target_global_axis_y[self.map_mpi[self.rank]["dst_global_y"]]
            else:
                return self.target_global_axis_y[self.map_mpi[self.rank]["dst_global_y"],self.map_mpi[self.rank]["dst_global_x"]]
        
    def find_point_index(self,target_lon, target_lat, decimal_tolerance=5, method="classic", only_mask_value=True,type="source"):
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
        lon = self.read_axis_x(type="source",with_overlap=False)
        lat = self.read_axis_y(type="source",with_overlap=False)

        if self.check_point_is_inside(target_lon, target_lat, lon, lat, tolerance=decimal_tolerance):
            try:
                mask = self.read_variable_2D_sea_binary_mask(type="source",with_overlap=False)
                print(np.shape(mask))
            except NotImplementedError:
                logging.warning("No 2D sea binary mask found")
                #mask = np.ones([self.source_global_y_size, self.source_global_x_size])
                only_mask_value = False

            dist = np.zeros([self.get_y_size(type="source",with_overlap=False), self.get_x_size(type="source",with_overlap=False)])
            dist[:] = 10000000

            if method=="classic":
                for x in range(0, self.get_x_size(type="source",with_overlap=False)):
                    for y in range(0, self.get_y_size(type="source",with_overlap=False)):

                        if only_mask_value:
                            if(mask[y,x] == 1): #=Terre
                                if self.reader.is_regular_grid():
                                    dist[y,x] = distance_on_unit_sphere(target_lon,target_lat,lon[x],lat[y])
                                else:
                                    dist[y,x] = distance_on_unit_sphere(target_lon,target_lat,lon[y,x],lat[y,x])
                        else:
                            if self.reader.is_regular_grid():
                                dist[y, x] = distance_on_unit_sphere(target_lon, target_lat, lon[x], lat[y])
                            else:
                                dist[y, x] = distance_on_unit_sphere(target_lon, target_lat, lon[y, x], lat[y, x])

                nearest_y_index,nearest_x_index = np.where(dist == np.min(dist))

                if len(nearest_y_index) == 0 or len(nearest_x_index) == 0:
                    logging.error("no neaarest point found")
                    raise RuntimeError("No nearest point found")

                nearest_x_index = nearest_x_index[0]
                nearest_y_index = nearest_y_index[0]
                min_dist = dist[nearest_y_index,nearest_x_index]

                if self.is_regular_grid(type="source"):
                    nearest_lon = lon[nearest_x_index]
                    nearest_lat = lat[nearest_y_index]
                else:
                    nearest_lon = lon[nearest_y_index, nearest_x_index]
                    nearest_lat = lat[nearest_y_index, nearest_x_index]

                if type == "source":
                    return [nearest_x_index, nearest_y_index, nearest_lon, nearest_lat, min_dist]
                elif type == "source_global":
                    return [self.map_mpi[self.rank]["src_global_x"].start+nearest_x_index, self.map_mpi[self.rank]["src_global_y"].start+nearest_y_index, nearest_lon, nearest_lat, min_dist]
                else:
                    raise ValueError("Type doesn't match [source, source_global]")

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

                    if type == "source":
                        return [nearest_x_index, nearest_y_index, nearest_lon, nearest_lat, min_dist]
                    elif type == "source_global":
                        return [self.map_mpi[self.rank]["src_global_x"].start + nearest_x_index,
                                self.map_mpi[self.rank]["src_global_y"].start + nearest_y_index, nearest_lon,
                                nearest_lat, min_dist]
                    else:
                        raise ValueError("Type doesn't match [source, source_global]")

                else:
                    raise NotImplementedError("Method " + str(method) + " is not implemented for regular grid.")

            else:
                raise RuntimeError("Method "+str(method)+" is not implemented yet.")
        else:
            logging.warning("Point is outside the rank n°"+str(self.rank))
            raise NotFoundInRankError(self.rank,"Point is outside the rank")
    
    # Variables
    #################
    # HYDRO
    # 2D
    #################
    def read_variable_bathymetry(self):     
        """Retourne la bathymétrie sur toute la couverture
    @return: un tableau en deux dimensions [y,x]."""
        data = self.reader.read_variable_bathymetry(
            self.map_mpi[self.rank]["src_global_x_overlap"].start,
            self.map_mpi[self.rank]["src_global_x_overlap"].stop,
            self.map_mpi[self.rank]["src_global_y_overlap"].start,
            self.map_mpi[self.rank]["src_global_y_overlap"].stop)

        if self.horizontal_resampling:
            data = resample_2d_to_grid(self.read_axis_x(type="source", with_overlap=True),
                                       self.read_axis_y(type="source", with_overlap=True),
                                       self.read_axis_x(type="target", with_overlap=True),
                                       self.read_axis_y(type="target", with_overlap=True),
                                       data,
                                       Coverage.HORIZONTAL_INTERPOLATION_METHOD)

        return data[self.map_mpi[self.rank]["dst_local_y"], self.map_mpi[self.rank]["dst_local_x"]]

    def read_variable_topography(self):
        """Retourne la topographie sur toute la couverture
    @return: un tableau en deux dimensions [y,x]."""
        data = self.reader.read_variable_topography(
            self.map_mpi[self.rank]["src_global_x_overlap"].start,
            self.map_mpi[self.rank]["src_global_x_overlap"].stop,
            self.map_mpi[self.rank]["src_global_y_overlap"].start,
            self.map_mpi[self.rank]["src_global_y_overlap"].stop)

        if self.horizontal_resampling:
            data = resample_2d_to_grid(self.read_axis_x(type="source", with_overlap=True),
                                       self.read_axis_y(type="source", with_overlap=True),
                                       self.read_axis_x(type="target", with_overlap=True),
                                       self.read_axis_y(type="target", with_overlap=True),
                                       data,
                                       Coverage.HORIZONTAL_INTERPOLATION_METHOD)

        return data[self.map_mpi[self.rank]["dst_local_y"], self.map_mpi[self.rank]["dst_local_x"]]
    
    def read_variable_mesh_size(self):     
        """Retourne la taille de la grille sur toute la couverture
    @return: un tableau en deux dimensions [y,x]."""
        data = self.reader.read_variable_mesh_size(
            self.map_mpi[self.rank]["src_global_x_overlap"].start,
            self.map_mpi[self.rank]["src_global_x_overlap"].stop,
            self.map_mpi[self.rank]["src_global_y_overlap"].start,
            self.map_mpi[self.rank]["src_global_y_overlap"].stop)

        if self.horizontal_resampling:
            data = resample_2d_to_grid(self.read_axis_x(type="source", with_overlap=True),
                                       self.read_axis_y(type="source", with_overlap=True),
                                       self.read_axis_x(type="target", with_overlap=True),
                                       self.read_axis_y(type="target", with_overlap=True),
                                       data,
                                       Coverage.HORIZONTAL_INTERPOLATION_METHOD)

        return data[self.map_mpi[self.rank]["dst_local_y"], self.map_mpi[self.rank]["dst_local_x"]]

    def read_variable_x_mesh_size(self):
        """Retourne la taille de la grille sur l'axe X
    @return: un tableau en deux dimensions [y,x]."""
        data = self.reader.read_variable_x_mesh_size(
            self.map_mpi[self.rank]["src_global_x_overlap"].start,
            self.map_mpi[self.rank]["src_global_x_overlap"].stop,
            self.map_mpi[self.rank]["src_global_y_overlap"].start,
            self.map_mpi[self.rank]["src_global_y_overlap"].stop)

        if self.horizontal_resampling:
            data = resample_2d_to_grid(self.read_axis_x(type="source", with_overlap=True),
                                       self.read_axis_y(type="source", with_overlap=True),
                                       self.read_axis_x(type="target", with_overlap=True),
                                       self.read_axis_y(type="target", with_overlap=True),
                                       data,
                                       Coverage.HORIZONTAL_INTERPOLATION_METHOD)

        return data[self.map_mpi[self.rank]["dst_local_y"], self.map_mpi[self.rank]["dst_local_x"]]

    def read_variable_y_mesh_size(self):
        """Retourne la taille de la grille sur l'axe Y
    @return: un tableau en deux dimensions [y,x]."""
        data = self.reader.read_variable_y_mesh_size(
            self.map_mpi[self.rank]["src_global_x_overlap"].start,
            self.map_mpi[self.rank]["src_global_x_overlap"].stop,
            self.map_mpi[self.rank]["src_global_y_overlap"].start,
            self.map_mpi[self.rank]["src_global_y_overlap"].stop)

        if self.horizontal_resampling:
            data = resample_2d_to_grid(self.read_axis_x(type="source", with_overlap=True),
                                       self.read_axis_y(type="source", with_overlap=True),
                                       self.read_axis_x(type="target", with_overlap=True),
                                       self.read_axis_y(type="target", with_overlap=True),
                                       data,
                                       Coverage.HORIZONTAL_INTERPOLATION_METHOD)

        return data[self.map_mpi[self.rank]["dst_local_y"], self.map_mpi[self.rank]["dst_local_x"]]
    
    def read_variable_2D_sea_binary_mask(self,type="target",with_overlap=False):
        """Retourne le masque terre/mer sur toute la couverture
    @return: un tableau en deux dimensions [y,x].
            0 = Terre
            1 = Mer
    """
        data = self.reader.read_variable_2D_sea_binary_mask(
            self.map_mpi[self.rank]["src_global_x_overlap"].start,
            self.map_mpi[self.rank]["src_global_x_overlap"].stop,
            self.map_mpi[self.rank]["src_global_y_overlap"].start,
            self.map_mpi[self.rank]["src_global_y_overlap"].stop)

        if type == "source" and with_overlap is True:
            return data[self.map_mpi[self.rank]["src_local_y_overlap"], self.map_mpi[self.rank]["src_local_x_overlap"]]
        elif type == "source" and with_overlap is False:
            return data[self.map_mpi[self.rank]["src_local_y"], self.map_mpi[self.rank]["src_local_x"]]
        elif type == "target":

            if self.horizontal_resampling:

                data = resample_2d_to_grid(self.read_axis_x(type="source", with_overlap=True),
                                           self.read_axis_y(type="source", with_overlap=True),
                                           self.read_axis_x(type="target", with_overlap=True),
                                           self.read_axis_y(type="target", with_overlap=True),
                                           data,
                                           Coverage.HORIZONTAL_INTERPOLATION_METHOD)

            return data[self.map_mpi[self.rank]["dst_local_y"],self.map_mpi[self.rank]["dst_local_x"]]

    def read_variable_Ha(self):
        """Retourne l'amplitude de la réanalyse
    @return: un tableau en deux dimensions [y,x]."""
        data = self.reader.read_variable_Ha(
            self.map_mpi[self.rank]["src_global_x_overlap"].start,
            self.map_mpi[self.rank]["src_global_x_overlap"].stop,
            self.map_mpi[self.rank]["src_global_y_overlap"].start,
            self.map_mpi[self.rank]["src_global_y_overlap"].stop)

        if self.horizontal_resampling:
            data = resample_2d_to_grid(self.read_axis_x(type="source", with_overlap=True),
                                       self.read_axis_y(type="source", with_overlap=True),
                                       self.read_axis_x(type="target", with_overlap=True),
                                       self.read_axis_y(type="target", with_overlap=True),
                                       data,
                                       Coverage.HORIZONTAL_INTERPOLATION_METHOD)

        return data[self.map_mpi[self.rank]["dst_local_y"], self.map_mpi[self.rank]["dst_local_x"]]




