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
from spatialetl.point.MultiPoint import MultiPoint
from spatialetl.operator.interpolator.InterpolatorCore import vertical_interpolation
import numpy as np
from datetime import datetime,timedelta
import logging

class LevelMultiPoint(MultiPoint):
    """"""
    DEPTH_DELTA = 1.0;  # meters

    #VERTICAL_INTERPOLATION_METHOD = None
    VERTICAL_INTERPOLATION_METHOD = "nearest"
    # VERTICAL_INTERPOLATION_METHOD = "pad"
    # VERTICAL_INTERPOLATION_METHOD = "backfill"
    # VERTICAL_INTERPOLATION_METHOD = "linear"

    def __init__(self,myReader,depths=None):
        MultiPoint.__init__(self, myReader)
        self.raw_levels = self.reader.read_axis_z()

        if self.is_sigma_coordinate(raw=True):
            self.raw_zmax = np.shape(self.raw_levels)[1]
        else:
            self.raw_zmax = np.shape(self.raw_levels)[0]

        if depths is None:
            self.target_levels =  self.raw_levels
            self.target_zmax = self.raw_zmax
        else:
            self.target_levels = np.array(depths)
            self.target_zmax = np.shape(self.target_levels)[0]

    # Axis
    def read_axis_z(self,raw=False):
        if raw:
            return self.raw_levels
        else:
            return self.target_levels

    def get_z_size(self,raw=False):
        if raw:
            return self.raw_zmax
        else:
            return self.target_zmax

    def is_sigma_coordinate(self,raw=False):

        if raw:
            if self.raw_levels.ndim == 2:
                return True
            elif self.raw_levels.ndim == 1:
                return False
            else:
                raise ValueError("Unable to recognize the type of the vertical grid.")
        else:
            if self.target_levels.ndim == 2:
                return True
            elif self.target_levels.ndim == 1:
                return False
            else:
                raise ValueError("Unable to recognize the type of the vertical grid.")


    def find_level_index(self, depth):
        """Retourne l'index de la profondeur la plus proche selon le point le plus proche.
    @type depth : integer ou flottant
    @param depth: Profondeur en mètre souhaitée ou index de la profondeur souhaitée
    @return:  un tableau de l'indice de la couche verticale inférieur la plus proche en chacun point de la grille. z < vert_coord[y,x] et z > vert_coord[y,x]+1.
    Les valeurs masquées valent -999."""

        vert_coord = np.empty([self.get_nb_points()], dtype=object)
        indexes_z = []

        if type(depth) == int or type(depth) == np.int32 or type(depth) == np.int64:

            if depth < 0 or depth >= self.get_z_size(raw=True):
                raise ValueError("Depth index have to range between 0 and " + str(
                    self.get_z_size(raw=True) - 1) + ". Actually depth index = " + str(depth))

            for i in range(0, self.get_nb_points()):
                    vert_coord[i] = []
                    vert_coord[i].append((depth))

            indexes_z.append((int(depth)))

        else:

            logging.debug("[LevelMultiPoint][find_level_index()] Looking for : " + str(depth)+" m water depth")

            if self.is_sigma_coordinate(raw=True) == True:  # Cas de grille sigma

                # Pour chaque point
                for i in range(0, self.get_nb_points()):

                    vert_coord[i] = []

                    logging.debug("[LevelMultiPoint][find_level_index()] Point "+str(i)+" - water depth candidates are : " + str(self.raw_levels[i,:]))

                    # On cherche l'index le plus proche
                    array = np.asarray(self.raw_levels[i,:])
                    index_z = (np.abs(array - depth)).argmin()

                    if abs(depth - self.raw_levels[i,index_z]) <= LevelMultiPoint.DEPTH_DELTA:
                        # On a trouvé une profondeur qui correspond au delta près.

                        logging.debug("[LevelMultiPoint][find_level_index()] Point "+str(i)+" - found : " + str(self.raw_levels[i,index_z])+" m water depth")

                        vert_coord[i].append((int(index_z)))
                        indexes_z.append((int(index_z)))

                        if abs(depth - self.raw_levels[i,index_z]) != 0.0:
                            # On n'a pas trouvé exactement notre profondeur, on va chercher autour au DEPTH_DELTA près
                            zz = index_z

                            while zz - 1 >= 0 and abs(self.raw_levels[i,zz - 1] - depth) <= LevelMultiPoint.DEPTH_DELTA and zz - 1 not in vert_coord[i]:
                                logging.debug("[LevelMultiPoint][find_level_index()] Point "+str(i)+" - found : " + str(self.raw_levels[i,zz-1])+" m water depth")
                                vert_coord[i].append((int(zz - 1)))
                                indexes_z.append((int(zz - 1)))
                                zz = zz - 1

                            zz = index_z
                            while zz + 1 < self.get_z_size(raw=True) and abs(self.raw_levels[i,zz + 1] - depth) <= LevelMultiPoint.DEPTH_DELTA and zz + 1 not in vert_coord[i]:
                                logging.debug("[LevelMultiPoint][find_level_index()] Point "+str(i)+" - found : " + str(self.raw_levels[i,zz+1])+" m water depth")
                                vert_coord[i].append((int(zz + 1)))
                                indexes_z.append((int(zz + 1)))
                                zz = zz + 1
                    else:
                        # raise ValueError("[LevelCoverage] " + str(depth) + " m water depth was not found for the point [" + str(x) + "," + str(y) + "]. Max depth found is for this point is " + str(self.levels[z,y,x]))
                        logging.debug(
                            "[LevelMultiPoint][find_level_index()] " + str(depth) + " m water depth was not found for the point [" + str(i) + "]. Max depth found is for this point is " + str(
                                self.raw_levels[i,index_z]) + " m.")

                    vert_coord[i] = np.array(np.unique(vert_coord[i]))

            else:  # Cas de grille classique

                logging.debug("[LevelMultiPoint][find_level_index()] Water depth candidates are : " + str(self.raw_levels))

                # On cherche l'index le plus proche
                array = np.asarray(self.raw_levels)
                index_z = (np.abs(array - depth)).argmin()

                if abs(depth - self.raw_levels[index_z]) <= LevelMultiPoint.DEPTH_DELTA:
                    # On a trouvé une profondeur qui correspond au delta près.

                    for i in range(0, self.get_nb_points()):

                        vert_coord[i] = []

                        logging.debug("[LevelMultiPoint][find_level_index()] Point "+str(i)+" - found : " + str(self.raw_levels[index_z]) + " m water depth")
                        vert_coord[i].append((int(index_z)))
                        indexes_z.append((int(index_z)))

                        if abs(depth - self.raw_levels[index_z]) != 0.0:
                            # On n'a pas trouvé exactement notre profondeur, on va chercher autour au DEPTH_DELTA près

                            zz = index_z
                            while zz - 1 >= 0 and abs(self.raw_levels[zz - 1] - depth) <= LevelMultiPoint.DEPTH_DELTA and zz - 1 not in vert_coord[i]:
                                logging.debug("[LevelMultiPoint][find_level_index()] Point "+str(i)+" - found : " + str(self.raw_levels[zz-1])+" m water depth")
                                vert_coord[i].append((int(zz - 1)))
                                indexes_z.append((int(zz - 1)))
                                zz = zz - 1

                            zz = index_z
                            while zz + 1 < self.get_z_size(raw=True) and abs(self.raw_levels[zz + 1] - depth) <= LevelMultiPoint.DEPTH_DELTA and zz + 1 not in vert_coord[i]:
                                logging.debug("[LevelMultiPoint][find_level_index()] Point "+str(i)+" - found : " + str(self.raw_levels[zz+1])+" m water depth")
                                vert_coord[i].append((int(zz + 1)))
                                indexes_z.append((int(zz + 1)))
                                zz = zz + 1

                else:
                    logging.debug("[LevelMultiPoint][find_level_index()] " + str(depth) + " m water depth was not found on the Z axis.")

        if not indexes_z:
            raise ValueError(
                "[LevelMultiPoint][find_level_index()] " + str(
                    depth) + " m water depth was not found in the grid. Maybe the LevelMultiPoint.DEPTH_DELTA (+/- " + str(
                    LevelMultiPoint.DEPTH_DELTA) + " m) is too small or the depth is out of range.")


        # On enlève les doublons
        result = np.array(np.unique(indexes_z))
        logging.debug("[LevelMultiPoint][find_level_index()] " + str(len(result)) + " vertical layer(s) found.")

        # On retourne le tableau d'index
        return [vert_coord, result]

    def interpolate_vertical(self,depth,vert_coord,indexes_z,layers):

        targetDepth = [depth]

        results = np.zeros([self.get_nb_points()])
        results[:] = np.NAN

        for i in range(0, self.get_nb_points()):

            if len(vert_coord[i]) == 1:
                # Il n'y a qu'une seule couche de sélectionner donc pas d'interpolation possible

                # On retrouve l'index de la layer
                array = np.asarray(indexes_z)
                index_layer = (np.abs(array - vert_coord[i][0])).argmin()

                results[i] = layers[index_layer, i]

            elif len(vert_coord[i]) > 1:

                candidateValues = np.zeros([len(vert_coord[i])])
                candidateDepths = np.zeros([len(vert_coord[i])])

                for z in range(0, len(vert_coord[i])):

                    # On retrouve l'index de la layer
                    array = np.asarray(indexes_z)
                    index_layer = (np.abs(array - vert_coord[i][z])).argmin()

                    if self.is_sigma_coordinate(raw=True):
                        candidateDepths[z] = self.raw_levels[i, indexes_z[index_layer]]
                    else:
                        candidateDepths[z] = self.raw_levels[indexes_z[index_layer]]

                    candidateValues[z] = layers[index_layer, i]

                results[i] = vertical_interpolation(candidateDepths, targetDepth, candidateValues,
                                                       LevelMultiPoint.VERTICAL_INTERPOLATION_METHOD)

        return results

    # Scalar
    def read_variable_sea_water_temperature_at_depth(self,depth):

        index_z = self.find_level_index(depth)

        #TODO
        # Il n'y a pas d'interpolation verticale pour le moment. Seule le module Coverage offre l'interpolation verticale.
        # On teste si le reader se base sur une Coverage, dans ce cas, on fournit directement la profondeur pour que le Coverage
        # interpolle verticalement pour nous.
        if self.is_coverage_based():
            index_z = depth

        data = self.reader.read_variable_sea_water_temperature_at_depth(index_z)

        return data

    def read_variable_sea_water_salinity_at_depth(self, depth):

        index_z = self.find_level_index(depth)

        #TODO
        # Il n'y a pas d'interpolation verticale pour le moment. Seule le module Coverage offre l'interpolation verticale.
        # On teste si le reader se base sur une Coverage, dans ce cas, on fournit directement la profondeur pour que le Coverage
        # interpolle verticalement pour nous.
        if self.is_coverage_based():
            index_z = depth

        data = self.reader.read_variable_sea_water_salinity_at_depth(index_z)

        return data

    def read_variable_sea_water_density_at_depth(self, depth):

        index_z = self.find_level_index(depth)

        #TODO
        # Il n'y a pas d'interpolation verticale pour le moment. Seule le module Coverage offre l'interpolation verticale.
        # On teste si le reader se base sur une Coverage, dans ce cas, on fournit directement la profondeur pour que le Coverage
        # interpolle verticalement pour nous.
        if self.is_coverage_based():
            index_z = depth

        data = self.reader.read_variable_sea_water_density_at_depth(index_z)

        return data

    def read_variable_sea_water_turbidity_at_depth(self, depth):

        index_z = self.find_level_index(depth)

        #TODO
        # Il n'y a pas d'interpolation verticale pour le moment. Seule le module Coverage offre l'interpolation verticale.
        # On teste si le reader se base sur une Coverage, dans ce cas, on fournit directement la profondeur pour que le Coverage
        # interpolle verticalement pour nous.
        if self.is_coverage_based():
            index_z = depth

        data = self.reader.read_variable_sea_water_turbidity_at_depth(index_z)

        return data

    def read_variable_sea_water_electrical_conductivity_at_depth(self, depth):

        index_z = self.find_level_index(depth)

        #TODO
        # Il n'y a pas d'interpolation verticale pour le moment. Seule le module Coverage offre l'interpolation verticale.
        # On teste si le reader se base sur une Coverage, dans ce cas, on fournit directement la profondeur pour que le Coverage
        # interpolle verticalement pour nous.
        if self.is_coverage_based():
            index_z = depth

        data = self.reader.read_variable_sea_water_electrical_conductivity_at_depth(index_z)

        return data


