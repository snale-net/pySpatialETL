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
from spatialetl.coverage.Coverage import Coverage
from scipy import spatial
import numpy as np
import logging

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
        self.target_global_z_size = self.source_global_z_size
        self.target_global_axis_z = self.source_global_axis_z

        if self.source_global_axis_z.ndim == 3:
            self.source_sigma_coordinate = True
        elif self.source_global_axis_z.ndim == 1:
            self.source_sigma_coordinate = False
        else:
            raise ValueError("[LevelCoverage] Unable to recognize the type of the vertical grid.")

        if resolution_z is not None:

            self.vertical_resampling = True
            self.target_sigma_coordinate = False

            if zbox == None:
                # we compute the destination grid
                zmin = np.min(self.source_global_axis_z)
                zmax = np.max(self.source_global_axis_z)
            else:
                zmin = zbox[0]
                zmax = zbox[1]

            self.target_global_axis_z = np.arange(zmin, zmax, resolution_z)
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

        else:
            self.target_global_axis_z = self.source_global_axis_z
            self.target_global_z_size = self.source_global_z_size
            self.target_sigma_coordinate = self.source_sigma_coordinate

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

        xmax=self.get_x_size(type="source",with_overlap=True)
        ymax=self.get_y_size(type="source",with_overlap=True)
        vert_coord = np.empty([ymax,xmax],dtype=object)
        indexes_z = []

        if depth in self.depth_weight:
            return self.depth_weight[depth]

        if type(depth) == int or type(depth) == np.int32 or type(depth) == np.int64:

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

                    if vert_coord[y,x] is None:
                        vert_coord[y, x] = []

                    vert_coord[y, x].append((int(index_z)))

                    if int(index_z) not in indexes_z:
                        indexes_z.append((int(index_z)))

            elif method=="classic":

                for y in range(0,ymax):
                    for x in range(0,xmax):

                        vert_coord[y, x] = []

                        if np.isnan(self.source_global_axis_z[:,y,x]).all() == False:

                            logging.debug("[LevelCoverage][find_level_index()] Point [" + str(x) + "," + str(y) + "] - water depth candidates are : " + str(self.source_global_axis_z[:,y,x]))

                            # On cherche l'index le plus proche
                            #array = np.asarray(self.source_global_axis_z[:,y,x])
                            #index_z = (np.abs(array - depth)).argmin()

                            if abs(depth - self.source_global_axis_z[index_z,y,x]) <= LevelCoverage.DEPTH_DELTA:
                                # On a trouvé une profondeur qui correspond au delta près.

                                logging.debug(
                                    "[LevelCoverage][find_level_index()] Point [" + str(x) + "," + str(y) + "] - found : " + str(
                                        self.source_global_axis_z[index_z, y, x]) + " m water depth")

                                vert_coord[y,x].append((int(index_z)))
                                indexes_z.append((int(index_z)))

                                if abs(depth - self.source_global_axis_z[index_z,y,x]) != 0.0:
                                    # On n'a pas trouvé exactement notre profondeur, on va chercher autour au DEPTH_DELTA près
                                    zz = index_z

                                    while zz - 1 >= 0 and abs(self.source_global_axis_z[zz-1,y,x]- depth) <= LevelCoverage.DEPTH_DELTA and zz -1 not in vert_coord[y,x]:
                                        logging.debug(
                                            "[LevelCoverage][find_level_index()] Point [" + str(x) + "," + str(y) + "] - found : " + str(
                                                self.source_global_axis_z[zz-1, y, x]) + " m water depth")
                                        vert_coord[y, x].append((int(zz-1)))
                                        indexes_z.append((int(zz-1)))
                                        zz = zz -1

                                    while zz + 1 < self.get_z_size() and abs(self.source_global_axis_z[zz+1,y,x]- depth) <= LevelCoverage.DEPTH_DELTA and zz +1 not in vert_coord[y,x]:
                                        logging.debug(
                                            "[LevelCoverage][find_level_index()] Point [" + str(x) + "," + str(y) + "] - found : " + str(
                                                self.source_global_axis_z[zz+1, y, x]) + " m water depth")
                                        vert_coord[y, x].append((int(zz+1)))
                                        indexes_z.append((int(zz+1)))
                                        zz = zz +1
                            else :
                                #raise ValueError("[LevelCoverage] " + str(depth) + " m water depth was not found for the point [" + str(x) + "," + str(y) + "]. Max depth found is for this point is " + str(self.levels[z,y,x]))
                                logging.debug("[LevelCoverage] " + str(depth) + " m water depth was not found for the point [" + str(x) + "," + str(y) + "]. Max depth found is for this point is " + str(self.source_global_axis_z[index_z,y,x])+" m.")

                            vert_coord[y,x] = np.array(np.unique(vert_coord[y,x]))

            else:
                raise ValueError("Unable to decode method between 'fast' or 'classic'.")

        else: # Cas de grille classique

            logging.debug("[LevelCoverage][find_level_index()] Looking for : " + str(depth) + " m water depth")

            # On cherche l'index le plus proche
            array = np.asarray(self.source_global_axis_z)
            index_z = (np.abs(array - depth)).argmin()

            if abs(depth - self.source_global_axis_z[index_z]) <= LevelCoverage.DEPTH_DELTA:
            # On a trouvé une profondeur qui correspond au delta près.

                logging.debug(
                "[LevelCoverage][find_level_index()] found : " + str(
                    self.source_global_axis_z[index_z]) + " m water depth")

                for y in range(0, ymax):
                    for x in range(0, xmax):

                        vert_coord[y, x] = []
                        vert_coord[y,x].append((int(index_z)))
                        indexes_z.append((int(index_z)))

                        if abs(depth - self.source_global_axis_z[index_z]) != 0.0:
                            # On n'a pas trouvé exactement notre profondeur, on va chercher autour au DEPTH_DELTA près
                            zz = index_z

                            while zz - 1 >= 0 and abs(self.source_global_axis_z[zz - 1] - depth) <= LevelCoverage.DEPTH_DELTA and zz -1 not in vert_coord[y,x]:
                                logging.debug(
                                    "[LevelCoverage][find_level_index()] found : " + str(
                                        self.source_global_axis_z[zz-1]) + " m water depth")

                                vert_coord[y,x].append((int(zz - 1)))
                                indexes_z.append((int(zz - 1)))
                                zz = zz - 1

                            while zz + 1 < self.get_z_size() and abs(self.source_global_axis_z[zz + 1] - depth) <= LevelCoverage.DEPTH_DELTA and zz +1 not in vert_coord[y,x]:
                                logging.debug(
                                    "[LevelCoverage][find_level_index()] found : " + str(
                                        self.source_global_axis_z[zz+1]) + " m water depth")
                                vert_coord[y,x].append((int(zz + 1)))
                                indexes_z.append((int(zz + 1)))
                                zz = zz + 1

            else:
                logging.debug("[LevelCoverage] " + str(depth) + " m water depth was not found on the Z axis.")

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
        
        
    

