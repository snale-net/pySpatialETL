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
from coverage.Coverage import Coverage
import numpy as np
import logging

class LevelCoverage(Coverage):    
    """
La classe LevelCoverage est une extension de la classe Coverage.
Elle rajoute une dimension verticale à la couverture horizontale classique.
"""
    DEPTH_DELTA = 0.5; #meters
    VERTICAL_INTERPOLATION_METHOD = "linear"
    
    def __init__(self, myReader):
        Coverage.__init__(self,myReader);
        self.levels = self.read_axis_z();
        self.last_index = None
        self.z_size = np.shape(self.levels)[0];
        self.depth_weight = []

        if self.levels.ndim == 3:
            self.sigma_coordinate = True
        elif self.levels.ndim == 1:
            self.sigma_coordinate = False
        else:
            raise ValueError("[LevelCoverage] Unable to recognize the type of the vertical grid.")
        
    # Axis
    def read_axis_z(self):
        """Retourne les valeurs (souvent en mètre) de l'axe z.
    @return:  si la grille est en coordonnée sigma alors un tableau à trois dimensions [z,y,x] est retourné sinon
    un tableau une dimension [z]."""
        return self.reader.read_axis_z();

    def is_sigma_coordinate(self):
        """Retourne vrai si la grille verticale est en coordonnée sigma, sinon faux.
    @return:  vrai si la grille verticale est en coordonnée sigma sinon faux."""
        return self.sigma_coordinate
    
    def get_z_size(self):
        """Retourne la taille de l'axe z.
    @return:  un entier correspondant à la taille de l'axe z."""
        return self.z_size
    
    def find_level_index(self,depth):
        """Retourne l'index de la profondeur la plus proche selon le point le plus proche.
    @type depth : integer ou flottant
    @param depth: Profondeur en mètre souhaitée ou index de la profondeur souhaitée
    @return:  un tableau de l'indice de la couche verticale inférieur la plus proche en chacun point de la grille. z < vert_coord[y,x] et z > vert_coord[y,x]+1.
    Les valeurs masquées valent -999."""

        xmax=self.get_x_size()
        ymax=self.get_y_size()
        zmax = self.get_z_size()
        #mask = self.read_variable_2D_sea_binary_mask()
        vert_coord = np.empty([ymax,xmax],dtype=object)
        indexes_z = []

        if type(depth) == int or type(depth) == np.int32 or type(depth) == np.int64:

            # A déplacer dans le reader
            #if depth < min or depth >= max:
            #    raise ValueError("Depth index have to range between "+str(min)+" and "+str(max)+". Actually depth index = "+str(depth))

            for y in range(0, ymax):
                for x in range(0, xmax):
                    vert_coord[y, x] = []
                    vert_coord[y,x].append((depth))

            indexes_z.append((int(depth)))

        elif self.is_sigma_coordinate() == True: # Cas de grille sigma

            logging.debug("[LevelCoverage][find_level_index()] Looking for : " + str(depth) + " m water depth")

            for y in range(0,ymax):
                for x in range(0,xmax):

                    vert_coord[y, x] = []

                    if np.isnan(self.levels[:,y,x]).all() == False:

                        logging.debug("[LevelCovergae][find_level_index()] Point [" + str(x) + "," + str(y) + "] - water depth candidates are : " + str(self.levels[:,y,x]))

                        # On cherche l'index le plus proche
                        array = np.asarray(self.levels[:,y,x])
                        index_z = (np.abs(array - depth)).argmin()

                        if abs(depth - self.levels[index_z,y,x]) <= LevelCoverage.DEPTH_DELTA:
                            # On a trouvé une profondeur qui correspond au delta près.

                            logging.debug(
                                "[LevelCoverage][find_level_index()] Point [" + str(x) + "," + str(y) + "] - found : " + str(
                                    self.levels[index_z, y, x]) + " m water depth")

                            vert_coord[y,x].append((int(index_z)))
                            indexes_z.append((int(index_z)))

                            if abs(depth - self.levels[index_z,y,x]) != 0.0:
                                # On n'a pas trouvé exactement notre profondeur, on va chercher autour au DEPTH_DELTA près
                                zz = index_z

                                while zz - 1 >= 0 and abs(self.levels[zz-1,y,x]- depth) <= LevelCoverage.DEPTH_DELTA and zz -1 not in vert_coord[y,x]:
                                    logging.debug(
                                        "[LevelCoverage][find_level_index()] Point [" + str(x) + "," + str(y) + "] - found : " + str(
                                            self.levels[zz-1, y, x]) + " m water depth")
                                    vert_coord[y, x].append((int(zz-1)))
                                    indexes_z.append((int(zz-1)))
                                    zz = zz -1

                                while zz + 1 < self.get_z_size() and abs(self.levels[zz+1,y,x]- depth) <= LevelCoverage.DEPTH_DELTA and zz +1 not in vert_coord[y,x]:
                                    logging.debug(
                                        "[LevelCoverage][find_level_index()] Point [" + str(x) + "," + str(y) + "] - found : " + str(
                                            self.levels[zz+1, y, x]) + " m water depth")
                                    vert_coord[y, x].append((int(zz+1)))
                                    indexes_z.append((int(zz+1)))
                                    zz = zz +1
                        else :
                            #raise ValueError("[LevelCoverage] " + str(depth) + " m water depth was not found for the point [" + str(x) + "," + str(y) + "]. Max depth found is for this point is " + str(self.levels[z,y,x]))
                            logging.debug("[LevelCoverage] " + str(depth) + " m water depth was not found for the point [" + str(x) + "," + str(y) + "]. Max depth found is for this point is " + str(self.levels[index_z,y,x])+" m.")

                        vert_coord[y,x] = np.array(np.unique(vert_coord[y,x]))

        else: # Cas de grille classique

            logging.debug("[LevelCoverage][find_level_index()] Looking for : " + str(depth) + " m water depth")

            # On cherche l'index le plus proche
            array = np.asarray(self.levels)
            index_z = (np.abs(array - depth)).argmin()

            if abs(depth - self.levels[index_z]) <= LevelCoverage.DEPTH_DELTA:
            # On a trouvé une profondeur qui correspond au delta près.

                logging.debug(
                "[LevelCoverage][find_level_index()] found : " + str(
                    self.levels[index_z]) + " m water depth")

                for y in range(0, ymax):
                    for x in range(0, xmax):

                        vert_coord[y, x] = []
                        vert_coord[y,x].append((int(index_z)))
                        indexes_z.append((int(index_z)))

                        if abs(depth - self.levels[index_z]) != 0.0:
                            # On n'a pas trouvé exactement notre profondeur, on va chercher autour au DEPTH_DELTA près
                            zz = index_z

                            while zz - 1 >= 0 and abs(self.levels[zz - 1] - depth) <= LevelCoverage.DEPTH_DELTA and zz -1 not in vert_coord[y,x]:
                                logging.debug(
                                    "[LevelCoverage][find_level_index()] found : " + str(
                                        self.levels[zz-1]) + " m water depth")

                                vert_coord[y,x].append((int(zz - 1)))
                                indexes_z.append((int(zz - 1)))
                                zz = zz - 1

                            while zz + 1 < self.get_z_size() and abs(self.levels[zz + 1] - depth) <= LevelCoverage.DEPTH_DELTA and zz +1 not in vert_coord[y,x]:
                                logging.debug(
                                    "[LevelCoverage][find_level_index()] found : " + str(
                                        self.levels[zz+1]) + " m water depth")
                                vert_coord[y,x].append((int(zz + 1)))
                                indexes_z.append((int(zz + 1)))
                                zz = zz + 1

            else:
                logging.debug("[LevelCoverage] " + str(depth) + " m water depth was not found on the Z axis.")

        if not indexes_z:
            raise ValueError(
                "[LevelCoverage] " + str(
                    depth) + " m water depth was not found in the grid. Maybe the LevelCoverage.DEPTH_DELTA (+/- " + str(
                    LevelCoverage.DEPTH_DELTA) + " m) is too small or the depth is out of range.")

        # On retourne le tableau d'index
        return [vert_coord,np.array(np.unique(indexes_z))]


    def read_variable_3D_sea_binary_mask(self):
        """Retourne le masque terre/mer sur toute la couverture selon la profondeur z
    @return: un tableau en deux dimensions [z,y,x].
            0 = Terre
            1 = Mer
    """
        return self.reader.read_variable_3D_sea_binary_mask()
        
        
    

