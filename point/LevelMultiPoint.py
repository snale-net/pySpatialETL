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
from point.MultiPoint import MultiPoint
import numpy as np
from datetime import datetime,timedelta
import logging

class LevelMultiPoint(MultiPoint):
    """"""
    DEPTH_DELTA = 1;  # meters

    def __init__(self,myReader):
        MultiPoint.__init__(self, myReader)
        self.levels = self.reader.read_axis_z()

    # Axis
    def read_axis_z(self):
        return self.reader.read_axis_z()

    def get_z_size(self):
        return np.shape(self.read_axis_z())[0]

    def find_level_index(self, depth):
        """Retourne l'index de la profondeur la plus proche selon le point le plus proche.
    @type depth : integer ou flottant
    @param depth: Profondeur en mètre souhaitée ou index de la profondeur souhaitée
    @return:  un tableau de l'indice de la couche verticale inférieur la plus proche en chacun point de la grille. z < vert_coord[y,x] et z > vert_coord[y,x]+1.
    Les valeurs masquées valent -999."""

        indexes_z = []

        if type(depth) == int:

            if depth < 0 or depth >= self.get_z_size():
                raise ValueError("Depth index have to range between 0 and " + str(
                    self.get_z_size() - 1) + ". Actually depth index = " + str(depth))

            indexes_z.append(int(depth));

        else:
            # On cherche l'index le plus proche
            array = np.asarray(self.levels)
            index_z = (np.abs(array - depth)).argmin()

            if abs(depth - self.levels[index_z]) <= LevelMultiPoint.DEPTH_DELTA:
                indexes_z.append(int(index_z))
            else:
                raise ValueError("[LevelMultiPoint] " + str(
                        depth) + " m water depth was not found. Maybe the LevelMultiPoint.DEPTH_DELTA (+/- " + str(
                        LevelMultiPoint.DEPTH_DELTA) + " m) is too small or the depth is out of range.")

        # On retourne le tableau d'index
        return np.array(np.unique(indexes_z))

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


