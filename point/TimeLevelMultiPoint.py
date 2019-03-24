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

from point.LevelMultiPoint import LevelMultiPoint
from point.TimeMultiPoint import TimeMultiPoint
from point.MultiPoint import MultiPoint
import numpy as np
from datetime import datetime,timedelta
import logging


class TimeLevelMultiPoint(LevelMultiPoint, TimeMultiPoint):
    """"""

    def __init__(self,myReader):
        LevelMultiPoint.__init__(self, myReader)
        TimeMultiPoint.__init__(self, myReader)

    # Scalar
    def read_variable_baroclinic_sea_water_velocity_at_time_and_depth(self,date,depth):

        index_t = self.find_time_index(date)
        index_z = self.find_level_index(depth)

        #TODO
        # Il n'y a pas d'interpolation verticale pour le moment. Seule le module Coverage offre l'interpolation verticale.
        # On teste si le reader se base sur une Coverage, dans ce cas, on fournit directement la profondeur pour que le Coverage
        # interpolle verticalement pour nous.
        if self.is_coverage_based():
            index_z = depth

        if len(index_t) > 1:
            layers = np.zeros([len(index_t), self.get_nb_points()])
            layers[::] = np.NAN

            for t in range(0, len(index_t)):
                layers[t] = self.reader.read_variable_baroclinic_sea_water_velocity_at_time_and_depth(index_t[t],index_z)

            data = self.interpolate_time(date, index_t, layers)

        else:
            data = self.reader.read_variable_baroclinic_sea_water_velocity_at_time_and_depth(index_t[0],index_z)

        return data

    def read_variable_sea_water_temperature_at_time_and_depth(self,date,depth):

        index_t = self.find_time_index(date)
        index_z = self.find_level_index(depth)

        #TODO
        # Il n'y a pas d'interpolation verticale pour le moment. Seule le module Coverage offre l'interpolation verticale.
        # On teste si le reader se base sur une Coverage, dans ce cas, on fournit directement la profondeur pour que le Coverage
        # interpolle verticalement pour nous.
        if self.is_coverage_based():
            index_z = depth

        if len(index_t) > 1:
            layers = np.zeros([len(index_t), self.get_nb_points()])
            layers[::] = np.NAN

            for t in range(0, len(index_t)):
                layers[t] = self.reader.read_variable_sea_water_temperature_at_time_and_depth(index_t[t],index_z)

            data = self.interpolate_time(date, index_t, layers)

        else:
            data = self.reader.read_variable_sea_water_temperature_at_time_and_depth(index_t[0],index_z)

        return data

    def read_variable_sea_water_salinity_at_time_and_depth(self, date, depth):

        index_t = self.find_time_index(date)
        index_z = self.find_level_index(depth)

        #TODO
        # Il n'y a pas d'interpolation verticale pour le moment. Seule le module Coverage offre l'interpolation verticale.
        # On teste si le reader se base sur une Coverage, dans ce cas, on fournit directement la profondeur pour que le Coverage
        # interpolle verticalement pour nous.
        if self.is_coverage_based():
            index_z = depth

        if len(index_t) > 1:
            layers = np.zeros([len(index_t), self.get_nb_points()])
            layers[::] = np.NAN

            for t in range(0, len(index_t)):
                layers[t] = self.reader.read_variable_sea_water_salinity_at_time_and_depth(index_t[t], index_z)

            data = self.interpolate_time(date, index_t, layers)

        else:
            data = self.reader.read_variable_sea_water_salinity_at_time_and_depth(index_t[0], index_z)

        return data


