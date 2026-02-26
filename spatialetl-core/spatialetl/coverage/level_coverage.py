#! /usr/bin/env python2.7
# -*- coding: utf-8 -*-
# MIT License
# Copyright (c) 2024 [SNALE - French SAS Company - RCS 951 724 616]
#
# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:
#
# The above copyright notice and this permission notice shall be included in all
# copies or substantial portions of the Software.
#
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
# SOFTWARE.
from __future__ import division, print_function, absolute_import

import inspect
import os
from itertools import product

import numpy as np

from spatialetl.coverage.coverage import Coverage
from spatialetl.operator.interpolator.interpolator_core import resample_2d_to_grid
from spatialetl.operator.interpolator.interpolator_core import vertical_interpolation
from spatialetl.utils.logger import logging


class LevelCoverage(Coverage):
    """
    La classe LevelCoverage est une extension de la classe Coverage.
    Elle rajoute une dimension verticale à la couverture horizontale classique.

    Examples
    --------
    >>> coverage = LevelCoverage(myReader, bbox=[-5, 5, -5, 5], resolution_x=0.1, resolution_y=0.1, resolution_z=0.1)
    """
    DEPTH_DELTA = 1.0;  # meters
    VERTICAL_INTERPOLATION_METHOD = "linear"

    def __init__(self, reader, bbox=None, resolution_x=None, resolution_y=None, zbox=None, resolution_z=None,
                 nb_thread: int = os.cpu_count() - 1):
        Coverage.__init__(self, reader, bbox=bbox, resolution_x=resolution_x, resolution_y=resolution_y,
                          nb_thread=nb_thread);

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
                raise ValueError("Zmin or Zmax is out of range. Current range: " + str(
                    np.min(self.source_global_axis_z)) + " / " + str(np.max(self.source_global_axis_z)))

            zmin = np.min(idx[0])
            zmax = np.max(idx[0]) + 1

            self.target_global_axis_z = self.source_global_axis_z[zmin:zmax]
            self.target_global_z_size = zmax - zmin
        else:
            idx = np.where((self.source_global_axis_z[:, 0:self.get_y_size(type="source_global"), 0:self.get_x_size(
                type="source_global")] >= Zmin) &
                           (self.source_global_axis_z[:, 0:self.get_y_size(type="source_global"), 0:self.get_x_size(
                               type="source_global")] <= Zmax))

            if (np.shape(idx)[1] == 0):
                raise ValueError("Zmin or Zmax is out of range. Current range: " + str(
                    np.min(self.source_global_axis_z)) + " / " + str(np.max(self.source_global_axis_z)))
            zmin = np.min(idx[0])
            zmax = np.max(idx[0]) + 1

            self.target_global_axis_z = self.source_global_axis_z[
                zmin:zmax, 0:self.get_y_size(type="source_global"), 0:self.get_x_size(type="source_global")]
            self.target_global_z_size = zmax - zmin

        # source_global sont réduit au zoom

        # On calcule la grille de destination
        self.target_sigma_coordinate = self.source_sigma_coordinate

        if resolution_z is not None:

            if resolution_z <= 0.0:
                raise ValueError("[LevelCoverage] resolution_z have to be upper than 0.0")

            self.vertical_resampling = True
            self.target_sigma_coordinate = False
            if resolution_z >= Zmax:  # On fait un seul niveau
                self.target_global_axis_z = [Zmax];
            else:
                self.target_global_axis_z = np.arange(Zmin, Zmax + resolution_z, resolution_z)

            self.target_global_z_size = len(self.target_global_axis_z)

            if self.rank == 0:

                if self.is_sigma_coordinate(type="source"):
                    logging.info(
                        '[vertical_interpolation] Source grid size : ' + str(
                            self.source_global_z_size) + " sigma coordinates level(s)")
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
    def read_axis_z(self, type="target", with_horizontal_overlap=False):
        """Retourne les valeurs (souvent en mètre) de l'axe z.
    @return:  si la grille est en coordonnée sigma alors un tableau à trois dimensions [z,y,x] est retourné sinon
    un tableau une dimension [z]."""

        if type == "source" and with_horizontal_overlap is True:

            if self.is_sigma_coordinate(type):
                return self.source_global_axis_z[:, self.parallel_map[self.rank]["src_global_y_overlap"],
                self.parallel_map[self.rank]["src_global_x_overlap"]]
            else:
                return self.source_global_axis_z

        elif type == "source" and with_horizontal_overlap is False:

            if self.is_sigma_coordinate(type):
                return self.source_global_axis_z[:, self.parallel_map[self.rank]["src_global_y"],
                self.parallel_map[self.rank]["src_global_x"]]
            else:
                return self.source_global_axis_z

        elif type == "target" and with_horizontal_overlap is True:

            if self.is_sigma_coordinate(type):
                return self.target_global_axis_z[:, self.parallel_map[self.rank]["dst_global_y_overlap"],
                self.parallel_map[self.rank]["dst_global_x_overlap"]]
            else:
                return self.target_global_axis_z
        else:

            if self.is_sigma_coordinate(type):
                return self.target_global_axis_z[:, self.parallel_map[self.rank]["dst_global_y"],
                self.parallel_map[self.rank]["dst_global_x"]]
            else:
                return self.target_global_axis_z

    def is_sigma_coordinate(self, type="target"):
        """Retourne vrai si la grille verticale est en coordonnée sigma, sinon faux.
    @return:  vrai si la grille verticale est en coordonnée sigma sinon faux."""
        if type == "source":
            return self.source_sigma_coordinate
        else:
            return self.target_sigma_coordinate

    def get_z_size(self, type="target"):
        """Retourne la taille de l'axe z.
    @return:  un entier correspondant à la taille de l'axe z."""

        if type == "source":
            return self.source_global_z_size
        else:
            return self.target_global_z_size

    def find_level_index(self, depth, method="fast"):
        """Retourne l'index de la profondeur la plus proche selon le point le plus proche.
    @type depth : integer ou flottant
    @param depth: Profondeur en mètre souhaitée ou index de la profondeur souhaitée
    @return:  un tableau de l'indice de la couche verticale inférieur la plus proche en chacun point de la grille. z < vert_coord[y,x] et z > vert_coord[y,x]+1.
    Les valeurs masquées valent -999."""

        if depth in self.depth_weight:
            return self.depth_weight[depth]

        xmax = self.get_x_size(type="source", with_overlap=True)
        ymax = self.get_y_size(type="source", with_overlap=True)
        vert_coord = np.empty([ymax, xmax], dtype=object)
        indexes_z = []

        if type(depth) == int or type(depth) == np.int32 or type(depth) == np.int64:

            if depth < 0 or depth >= self.get_z_size(type="source"):
                raise ValueError("Depth index have to range between 0 and " + str(
                    self.get_z_size(type="source") - 1) + ". Actually Depth index = " + str(depth))

            for y in range(0, ymax):
                for x in range(0, xmax):
                    vert_coord[y, x] = []
                    vert_coord[y, x].append((depth))

            indexes_z.append((int(depth)))

        elif self.is_sigma_coordinate(type="source") == True:  # Cas de grille sigma

            if self.rank == 0:
                logging.debug("[LevelCoverage][find_level_index()] Looking for : " + str(
                    depth) + " m water depth with an interval of +/- " + str(LevelCoverage.DEPTH_DELTA) + " m")

            if method == "fast":

                X = np.abs(self.source_global_axis_z[
                               :, self.parallel_map[self.rank]["src_global_y_overlap"], self.parallel_map[self.rank][
                                   "src_global_x_overlap"]] - depth)
                idx = np.where(X <= LevelCoverage.DEPTH_DELTA)
                vert_coord[:] = None

                for index in range(np.shape(idx)[1]):
                    index_z = idx[0][index]
                    x = idx[2][index]
                    y = idx[1][index]

                    if vert_coord[y, x] is None:  # first time
                        vert_coord[y, x] = []

                    vert_coord[y, x].append((int(index_z)))

                    if int(index_z) not in indexes_z:
                        indexes_z.append((int(index_z)))
            else:
                raise ValueError("Unable to decode method between 'fast' or 'classic'.")

        else:  # Cas de grille classique

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

                    for y, x in product(range(0, ymax), range(0, xmax)):
                        if vert_coord[y, x] is None:  # first time
                            vert_coord[y, x] = []

                        vert_coord[y, x].append((int(index_z)))

                        if int(index_z) not in indexes_z:
                            indexes_z.append((int(index_z)))
                else:
                    idx = np.where(X <= LevelCoverage.DEPTH_DELTA)

                    for index in range(np.shape(idx)[1]):
                        index_z = idx[0][index]

                        for y, x in product(range(0, ymax), range(0, xmax)):
                            if vert_coord[y, x] is None:  # first time
                                vert_coord[y, x] = []

                            # logging.debug(
                            #    "[LevelCoverage][find_level_index()] found : " + str(
                            #        self.source_global_axis_z[index_z]) + " m water depth")

                            vert_coord[y, x].append((int(index_z)))

                            if int(index_z) not in indexes_z:
                                indexes_z.append((int(index_z)))

            else:
                raise ValueError("Unable to decode method between 'fast' or 'classic'.")

        if len(indexes_z) == 0:
            logging.warning("[LevelCoverage] " + str(
                depth) + " m water depth was not found in the grid (proc n° " + str(
                self.rank) + "). Maybe the LevelCoverage.DEPTH_DELTA (+/- " + str(
                LevelCoverage.DEPTH_DELTA) + " m) is too small or the depth is out of range.")

        if self.rank == 0:
            logging.debug("[LevelCoverage][find_level_index()] Found " + str(len(indexes_z)) + " candidate level(s)")

        # On retourne le tableau d'index
        self.depth_weight[depth] = [vert_coord, np.array(np.unique(indexes_z))]

        return self.depth_weight[depth]

    def __read_variable(self, function_name, depth):

        fn = getattr(self.reader, function_name)
        vert_coord, indexes_z = self.find_level_index(depth);
        self.layers_temp[::] = np.nan
        self.data_temp[::] = np.nan
        targetDepth = [depth]

        for z in range(0, len(indexes_z)):
            self.layers_temp[z] = fn(
                indexes_z[z],
                self.parallel_map[self.rank]["src_global_x_overlap"].start,
                self.parallel_map[self.rank]["src_global_x_overlap"].stop,
                self.parallel_map[self.rank]["src_global_y_overlap"].start,
                self.parallel_map[self.rank]["src_global_y_overlap"].stop)

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

                    candidateValues[z] = self.layers_temp[index_layer, 0, y, x]

                self.data_temp[0, y, x] = vertical_interpolation(candidateDepths, targetDepth, candidateValues,
                                                                 LevelCoverage.VERTICAL_INTERPOLATION_METHOD)

        is_vector = True if len(np.shape(self.data_temp[0])) == 3 and np.shape(self.data_temp[0])[0] == 2 else False

        if self.horizontal_resampling:
            if is_vector:
                return [self.resample_2d_variable(self.data_temp[0, 0]),
                        self.resample_2d_variable(self.data_temp[0, 1])]
            else:
                return self.resample_2d_variable(self.data_temp[0])
        else:
            if is_vector:
                return [
                    self.data_temp[0, 0, self.parallel_map[self.rank]["dst_global_y"], self.parallel_map[self.rank][
                        "dst_global_x"]],
                    self.data_temp[0, 1, self.parallel_map[self.rank]["dst_global_y"], self.parallel_map[self.rank][
                        "dst_global_x"]]
                ]
            else:
                return self.data_temp[
                    0, self.parallel_map[self.rank]["dst_global_y"], self.parallel_map[self.rank]["dst_global_x"]]

    def read_variable_3D_sea_binary_mask_at_depth(self, depth):
        """Retourne le masque terre/mer sur toute la couverture selon la profondeur z
    @return: un tableau en deux dimensions [z,y,x].
            0 = Terre
            1 = Mer
    """
        return self.__read_variable(inspect.stack()[0][3], depth=depth)

    def read_variable_depth_at_depth(self, depth):

        return self.__read_variable(inspect.stack()[0][3], depth=depth)
