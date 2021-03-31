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

from itertools import product

import numpy as np

from spatialetl.coverage.Coverage import Coverage
from spatialetl.operator.interpolator.InterpolatorCore import resample_2d_to_grid
from spatialetl.operator.interpolator.InterpolatorCore import vertical_interpolation
from spatialetl.utils.logger import logging


class LevelCoverage(Coverage):
    """
La classe LevelCoverage est une extension de la classe Coverage.
Elle rajoute une dimension verticale à la couverture horizontale classique.
"""
    DEPTH_DELTA = 1.0; #meters
    VERTICAL_INTERPOLATION_METHOD = "linear"
    
    def __init__(self, myReader,bbox=None,resolution_x=None,resolution_y=None,zbox=None,resolution_z=None):
        Coverage.__init__(self,myReader,bbox=bbox, resolution_x=resolution_x, resolution_y=resolution_y);

        self.vertical_resampling = False
        self.source_sigma_coordinate = False
        self.target_sigma_coordinate = False

        self.source_global_z_size = self.reader.get_z_size()
        self.source_global_axis_z = self.reader.read_axis_z();

        if self.source_global_axis_z.ndim == 3:
            self.source_sigma_coordinate = True
        elif self.source_global_axis_z.ndim == 1:
            self.source_sigma_coordinate = False
        else:
            raise ValueError("[LevelCoverage] Unable to recognize the type of the vertical grid.")

        if zbox == None:
            # we compute the destination grid
            Zmin = np.min(self.source_global_axis_z)
            Zmax = np.max(self.source_global_axis_z)
        else:
            Zmin = zbox[0]
            Zmax = zbox[1]

        if not self.source_sigma_coordinate:
            idx = np.where((self.source_global_axis_z >= Zmin) &
                           (self.source_global_axis_z <= Zmax))

            if (np.shape(idx)[1] == 0):
                raise ValueError("Zmin & Zmax out of range")

            zmin = np.min(idx[0])
            zmax = np.max(idx[0]) + 1

            self.target_global_axis_z = self.source_global_axis_z[zmin:zmax]
            self.target_global_z_size = zmax - zmin
        else:
            idx = np.where((self.source_global_axis_z[:,0:self.get_y_size(type="source_global"),0:self.get_x_size(type="source_global")] >= Zmin) &
                           (self.source_global_axis_z[:,0:self.get_y_size(type="source_global"),0:self.get_x_size(type="source_global")] <= Zmax))

            if (np.shape(idx)[1] == 0):
                raise ValueError("Zmin & Zmax out of range")
            zmin = np.min(idx[0])
            zmax = np.max(idx[0]) + 1

            self.target_global_axis_z = self.source_global_axis_z[zmin:zmax,0:self.get_y_size(type="source_global"),0:self.get_x_size(type="source_global")]
            self.target_global_z_size = zmax - zmin

        # source_global sont réduit au zoom

        # On calcule la grille de destination
        self.target_sigma_coordinate = self.source_sigma_coordinate

        if resolution_z is not None:

            if resolution_z <= 0.0:
                raise ValueError("[LevelCoverage] resolution_z have to be upper than 0.0")

            self.vertical_resampling = True
            self.target_sigma_coordinate = False
            self.target_global_axis_z = np.arange(Zmin, Zmax+resolution_z, resolution_z)
            self.target_global_z_size = len(self.target_global_axis_z)

            if self.rank == 0:

                if self.is_sigma_coordinate(type="source"):
                    logging.info(
                        '[vertical_interpolation] Source grid size : ' + str(self.source_global_z_size) + " sigma coordinates level(s)")
                else:
                    logging.info(
                        '[vertical_interpolation] Source grid size : ' + str(self.source_global_z_size) + " level(s)")

                logging.info(
                    '[vertical_interpolation] Target grid size : ' + str(self.target_global_z_size) + " level(s)")

        self.depth_weight = {}

        if type(self) == LevelCoverage and self.horizontal_resampling and self.rank == 0:
            logging.info(
                '[horizontal_interpolation] Source grid size : (' + str(self.source_global_x_size) + ", " + str(
                    self.source_global_y_size) + ")")
            logging.info(
                '[horizontal_interpolation] Target grid size : (' + str(self.target_global_x_size) + ", " + str(
                    self.target_global_y_size) + ")")

    # Axis
    def read_axis_z(self,type="target",with_horizontal_overlap=False):
        """Retourne les valeurs (souvent en mètre) de l'axe z.
    @return:  si la grille est en coordonnée sigma alors un tableau à trois dimensions [z,y,x] est retourné sinon
    un tableau une dimension [z]."""

        if type == "source" and with_horizontal_overlap is True:

            if self.is_sigma_coordinate(type):
                return self.source_global_axis_z[:,self.map_mpi[self.rank]["src_global_y_overlap"],
                                                 self.map_mpi[self.rank]["src_global_x_overlap"]]
            else:
                return self.source_global_axis_z

        elif type == "source" and with_horizontal_overlap is False:

            if self.is_sigma_coordinate(type):
                return self.source_global_axis_z[:,self.map_mpi[self.rank]["src_global_y"],
                                                 self.map_mpi[self.rank]["src_global_x"]]
            else:
                return self.source_global_axis_z

        elif type == "target" and with_horizontal_overlap is True:

            if self.is_sigma_coordinate(type):
                return self.target_global_axis_z[:,self.map_mpi[self.rank]["dst_global_y_overlap"],
                                                 self.map_mpi[self.rank]["dst_global_x_overlap"]]
            else:
                return self.target_global_axis_z
        else:

            if self.is_sigma_coordinate(type):
                return self.target_global_axis_z[:,self.map_mpi[self.rank]["dst_global_y"],
                                                 self.map_mpi[self.rank]["dst_global_x"]]
            else:
                return self.target_global_axis_z

    def is_sigma_coordinate(self,type="target"):
        """Retourne vrai si la grille verticale est en coordonnée sigma, sinon faux.
    @return:  vrai si la grille verticale est en coordonnée sigma sinon faux."""
        if type == "source":
            return self.source_sigma_coordinate
        else:
            return self.target_sigma_coordinate
    
    def get_z_size(self,type="target"):
        """Retourne la taille de l'axe z.
    @return:  un entier correspondant à la taille de l'axe z."""

        if type=="source":
            return self.source_global_z_size
        else:
            return self.target_global_z_size

    def find_level_index(self,depth,method="fast"):
        """Retourne l'index de la profondeur la plus proche selon le point le plus proche.
    @type depth : integer ou flottant
    @param depth: Profondeur en mètre souhaitée ou index de la profondeur souhaitée
    @return:  un tableau de l'indice de la couche verticale inférieur la plus proche en chacun point de la grille. z < vert_coord[y,x] et z > vert_coord[y,x]+1.
    Les valeurs masquées valent -999."""

        if depth in self.depth_weight:
            return self.depth_weight[depth]

        xmax=self.get_x_size(type="source",with_overlap=True)
        ymax=self.get_y_size(type="source",with_overlap=True)
        vert_coord = np.empty([ymax,xmax],dtype=object)
        indexes_z = []

        if type(depth) == int or type(depth) == np.int32 or type(depth) == np.int64:

            if depth < 0 or depth >= self.get_z_size(type="source"):
                raise ValueError("Depth index have to range between 0 and " + str(
                    self.get_z_size(type="source") - 1) + ". Actually Depth index = " + str(depth))

            for y in range(0, ymax):
                for x in range(0, xmax):
                    vert_coord[y, x] = []
                    vert_coord[y,x].append((depth))

            indexes_z.append((int(depth)))

        elif self.is_sigma_coordinate(type="source") == True: # Cas de grille sigma

            if self.rank == 0:
                logging.debug("[LevelCoverage][find_level_index()] Looking for : " + str(depth) + " m water depth with an interval of +/- " + str(LevelCoverage.DEPTH_DELTA) + " m")

            if method == "fast":

                X = np.abs(self.source_global_axis_z[:,self.map_mpi[self.rank]["src_global_y_overlap"],self.map_mpi[self.rank]["src_global_x_overlap"]] - depth)
                idx = np.where(X <= LevelCoverage.DEPTH_DELTA)
                vert_coord[:] = None

                for index in range(np.shape(idx)[1]):
                    index_z = idx[0][index]
                    x =  idx[2][index]
                    y =  idx[1][index]

                    if vert_coord[y,x] is None: #first time
                        vert_coord[y, x] = []

                    vert_coord[y, x].append((int(index_z)))

                    if int(index_z) not in indexes_z:
                        indexes_z.append((int(index_z)))
            else:
                raise ValueError("Unable to decode method between 'fast' or 'classic'.")

        else: # Cas de grille classique

            if self.rank == 0:
                logging.debug("[LevelCoverage][find_level_index()] Looking for : " + str(
                    depth) + " m water depth with an interval of +/- " + str(LevelCoverage.DEPTH_DELTA) + " m")

            if method == "fast":

                X = np.abs(self.source_global_axis_z - depth)

                idx = np.where(X == 0.0)
                if (len(idx[0]) == 1):
                    index_z = idx[0][0]

                    if self.rank == 0:
                        logging.debug("[LevelCoverage][find_level_index()] found : " + str(
                                self.source_global_axis_z[index_z]) + " m water depth")

                    for y,x in product(range(0,ymax),range(0,xmax)):
                            if vert_coord[y, x] is None:  # first time
                                vert_coord[y, x] = []

                            vert_coord[y, x].append((int(index_z)))

                            if int(index_z) not in indexes_z:
                                indexes_z.append((int(index_z)))
                else:
                    idx = np.where(X <= LevelCoverage.DEPTH_DELTA)

                    for index in range(np.shape(idx)[1]):
                        index_z = idx[0][index]

                        for y,x in product(range(0,ymax),range(0,xmax)):
                            if vert_coord[y, x] is None:  # first time
                                vert_coord[y, x] = []

                            #logging.debug(
                            #    "[LevelCoverage][find_level_index()] found : " + str(
                            #        self.source_global_axis_z[index_z]) + " m water depth")

                            vert_coord[y, x].append((int(index_z)))

                            if int(index_z) not in indexes_z:
                                indexes_z.append((int(index_z)))

            else:
                raise ValueError("Unable to decode method between 'fast' or 'classic'.")

        if len(indexes_z)==0:
            logging.warning("[LevelCoverage] " + str(
                    depth) + " m water depth was not found in the grid (proc n° "+str(self.rank)+"). Maybe the LevelCoverage.DEPTH_DELTA (+/- " + str(
                    LevelCoverage.DEPTH_DELTA) + " m) is too small or the depth is out of range.")

        if self.rank == 0:
            logging.debug("[LevelCoverage][find_level_index()] Found " + str(len(indexes_z)) + " candidate level(s)")

        # On retourne le tableau d'index
        self.depth_weight[depth] = [vert_coord,np.array(np.unique(indexes_z))]

        return self.depth_weight[depth]

    def read_variable_3D_sea_binary_mask(self):
        """Retourne le masque terre/mer sur toute la couverture selon la profondeur z
    @return: un tableau en deux dimensions [z,y,x].
            0 = Terre
            1 = Mer
    """
        return self.reader.read_variable_3D_sea_binary_mask()

    def read_variable_depth_at_depth(self, depth):
        vert_coord,indexes_z = self.find_level_index(depth);
        self.layers_temp[::] = np.NAN
        self.data_temp[::] = np.NAN
        targetDepth = [depth]

        for z in range(0, len(indexes_z)):
            self.layers_temp[z] = self.reader.read_variable_depth_at_depth(
                indexes_z[z],
                self.map_mpi[self.rank]["src_global_x_overlap"].start,
                self.map_mpi[self.rank]["src_global_x_overlap"].stop,
                self.map_mpi[self.rank]["src_global_y_overlap"].start,
                self.map_mpi[self.rank]["src_global_y_overlap"].stop)

        idx = np.where(vert_coord != None)
        for index in range(np.shape(idx)[1]):
            x = idx[1][index]
            y = idx[0][index]

            if len(vert_coord[y, x]) == 1:
                # Il n'y a qu'une seule couche de sélectionner donc pas d'interpolation possible
                # On retrouve l'index de la layer
                index_layer = (np.abs(indexes_z - vert_coord[y, x][0])).argmin()
                self.data_temp[0, y, x] = self.layers_temp[index_layer, 0, y, x]

            else:
                candidateValues = np.zeros([len(vert_coord[y, x])])
                candidateDepths = np.zeros([len(vert_coord[y, x])])

                for z in range(0, len(vert_coord[y, x])):
                    # On retrouve l'index de la layer
                    index_layer = (np.abs(indexes_z - vert_coord[y, x][z])).argmin()

                    if self.is_sigma_coordinate(type="source"):
                        candidateDepths[z] = self.read_axis_z(type="source", with_horizontal_overlap=True)[
                            vert_coord[y, x][z], y, x]
                    else:
                        candidateDepths[z] = self.read_axis_z(type="source", with_horizontal_overlap=True)[
                            vert_coord[y, x][z]]

                    candidateValues[z] = self.layers_temp[index_layer,0, y, x]

                self.data_temp[0,y, x] = vertical_interpolation(candidateDepths, targetDepth, candidateValues,
                                                    LevelCoverage.VERTICAL_INTERPOLATION_METHOD)

        if self.horizontal_resampling:
            return resample_2d_to_grid(self.read_axis_x(type="source", with_overlap=True),
                                       self.read_axis_y(type="source", with_overlap=True),
                                       self.read_axis_x(type="target", with_overlap=True),
                                       self.read_axis_y(type="target", with_overlap=True),
                                       self.data_temp[0],
                                       Coverage.HORIZONTAL_INTERPOLATION_METHOD)[self.map_mpi[self.rank]["dst_local_y"], self.map_mpi[self.rank]["dst_local_x"]]

        return self.data_temp[0,self.map_mpi[self.rank]["dst_local_y"], self.map_mpi[self.rank]["dst_local_x"]]
        
        
    

