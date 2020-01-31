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
from spatialetl.coverage.Coverage import  Coverage
from spatialetl.coverage.LevelCoverage import LevelCoverage
from spatialetl.coverage.TimeCoverage import TimeCoverage
from spatialetl.operator.interpolator.InterpolatorCore import resample_2d_to_grid
from spatialetl.operator.interpolator.InterpolatorCore import vertical_interpolation
import numpy as np
import logging

class TimeLevelCoverage(LevelCoverage,TimeCoverage):
    """La classe TimeLevelCoverage est une extension de la classe Coverage, LevelCoverage, TimeCoverage.
Elle rajoute les dimensions temporelle et verticale à la couverture horizontale classique.
    """
    def __init__(self, myReader,bbox=None,resolution_x=None,resolution_y=None,zbox=None,resolution_z=None,start_time=None,end_time=None,freq=None):

        LevelCoverage.__init__(self,myReader,bbox=bbox,resolution_x=resolution_x,resolution_y=resolution_y,zbox=zbox,resolution_z=resolution_z);
        TimeCoverage.__init__(self,myReader,bbox=bbox,resolution_x=resolution_x,resolution_y=resolution_y,start_time=start_time,end_time=end_time,freq=freq);

        if self.horizontal_resampling and self.rank == 0:
            logging.info(
                '[horizontal_interpolation] Source grid size : (' + str(self.source_global_x_size) + ", " + str(
                    self.source_global_y_size) + ")")
            logging.info(
                '[horizontal_interpolation] Target grid size : (' + str(self.target_global_x_size) + ", " + str(
                    self.target_global_y_size) + ")")

    #################
    # HYDRO
    # 3D
    #################

    def read_variable_sea_water_temperature_at_time_and_depth(self, time, depth):
        """Retourne la salinité à la date souhaitée et au niveau souhaité sur toute la couverture horizontale.
    @type time: datetime ou l'index
    @param time: date souhaitée
    @type depth: profondeur en mètre (float) ou index (integer)
    @param depth: profondeur souhaitée. Si le z est un entier, on considère qu'il s'agit de l'index,
    si c'est un flottant on considère qu'il s'agit d'une profondeur
    @return: un tableau en deux dimensions [y,x]."""

        index_t = self.find_time_index(time);
        tmp = self.find_level_index(depth);
        vert_coord = tmp[0]
        indexes_z = tmp[1]

        xmax = self.get_x_size(type="source", with_overlap=True)
        ymax = self.get_y_size(type="source", with_overlap=True)
        z_axis = self.read_axis_z(type="source", with_horizontal_overlap=True)

        layers = np.zeros([np.shape(indexes_z)[0], ymax, xmax])
        layers[::] = np.NAN

        data = np.zeros([ymax, xmax])
        data[:] = np.NAN
        targetDepth = [depth]

        mask_t = self.reader.read_variable_2D_sea_binary_mask(
                self.map_mpi[self.rank]["src_global_x_overlap"].start,
                self.map_mpi[self.rank]["src_global_x_overlap"].stop,
                self.map_mpi[self.rank]["src_global_y_overlap"].start,
                self.map_mpi[self.rank]["src_global_y_overlap"].stop)

        for z in range(0, len(indexes_z)):
            layers[z] = self.reader.read_variable_sea_water_temperature_at_time_and_depth(
                self.map_mpi[self.rank]["src_global_t"].start + index_t, indexes_z[z],
                self.map_mpi[self.rank]["src_global_x_overlap"].start,
                self.map_mpi[self.rank]["src_global_x_overlap"].stop,
                self.map_mpi[self.rank]["src_global_y_overlap"].start,
                self.map_mpi[self.rank]["src_global_y_overlap"].stop)

        idx = np.where(vert_coord != None)
        for index in range(np.shape(idx)[1]):
            x = idx[1][index]
            y = idx[0][index]

            if mask_t[y,x] == 1:

                if len(vert_coord[y, x]) == 1:
                    # Il n'y a qu'une seule couche de sélectionner donc pas d'interpolation possible
                    # On retrouve l'index de la layer
                    array = np.asarray(indexes_z)
                    index_layer = (np.abs(array - vert_coord[y, x][0])).argmin()
                    data[y, x] = layers[index_layer, y, x]

                elif len(vert_coord[y, x]) > 1:

                    candidateValues = np.zeros([len(vert_coord[y, x])])
                    candidateDepths = np.zeros([len(vert_coord[y, x])])

                    for z in range(0, len(vert_coord[y, x])):
                        # On retrouve l'index de la layer
                        array = np.asarray(indexes_z)
                        index_layer = (np.abs(array - vert_coord[y, x][z])).argmin()

                        if self.is_sigma_coordinate(type="source"):
                            candidateDepths[z] = z_axis[index_layer, y, x]
                        else:
                            candidateDepths[z] = z_axis[index_layer]

                        candidateValues[z] = layers[index_layer, y, x]

                    data[y, x] = vertical_interpolation(candidateDepths, targetDepth, candidateValues,
                                                        LevelCoverage.VERTICAL_INTERPOLATION_METHOD)

        if self.horizontal_resampling:
            data = resample_2d_to_grid(self.read_axis_x(type="source_global", with_overlap=True),
                                       self.read_axis_y(type="source_global", with_overlap=True),
                                       self.read_axis_x(type="target", with_overlap=True),
                                       self.read_axis_y(type="target", with_overlap=True),
                                       data,
                                       Coverage.HORIZONTAL_INTERPOLATION_METHOD)

        return data[self.map_mpi[self.rank]["dst_local_y"], self.map_mpi[self.rank]["dst_local_x"]]

    def read_variable_sea_water_salinity_at_time_and_depth(self, time, depth):
        """Retourne la salinité à la date souhaitée et au niveau souhaité sur toute la couverture horizontale.
    @type time: datetime ou l'index
    @param time: date souhaitée
    @type depth: profondeur en mètre (float) ou index (integer)
    @param depth: profondeur souhaitée. Si le z est un entier, on considère qu'il s'agit de l'index,
    si c'est un flottant on considère qu'il s'agit d'une profondeur
    @return: un tableau en deux dimensions [y,x]."""

        index_t = self.find_time_index(time);
        tmp = self.find_level_index(depth);
        vert_coord = tmp[0]
        indexes_z = tmp[1]

        xmax = self.get_x_size(type="source", with_overlap=True)
        ymax = self.get_y_size(type="source", with_overlap=True)
        z_axis = self.read_axis_z(type="source", with_horizontal_overlap=True)

        layers = np.zeros([np.shape(indexes_z)[0], ymax, xmax])
        layers[::] = np.NAN

        data = np.zeros([ymax, xmax])
        data[:] = np.NAN
        targetDepth = [depth]

        mask_t = self.reader.read_variable_2D_sea_binary_mask(
            self.map_mpi[self.rank]["src_global_x_overlap"].start,
            self.map_mpi[self.rank]["src_global_x_overlap"].stop,
            self.map_mpi[self.rank]["src_global_y_overlap"].start,
            self.map_mpi[self.rank]["src_global_y_overlap"].stop)

        for z in range(0, len(indexes_z)):
            layers[z] = self.reader.read_variable_sea_water_salinity_at_time_and_depth(
                self.map_mpi[self.rank]["src_global_t"].start + index_t, indexes_z[z],
                self.map_mpi[self.rank]["src_global_x_overlap"].start,
                self.map_mpi[self.rank]["src_global_x_overlap"].stop,
                self.map_mpi[self.rank]["src_global_y_overlap"].start,
                self.map_mpi[self.rank]["src_global_y_overlap"].stop)

        idx = np.where(vert_coord != None)
        for index in range(np.shape(idx)[1]):
            x = idx[1][index]
            y = idx[0][index]

            if mask_t[y, x] == 1:

                if len(vert_coord[y, x]) == 1:
                    # Il n'y a qu'une seule couche de sélectionner donc pas d'interpolation possible
                    # On retrouve l'index de la layer
                    array = np.asarray(indexes_z)
                    index_layer = (np.abs(array - vert_coord[y, x][0])).argmin()
                    data[y, x] = layers[index_layer, y, x]

                elif len(vert_coord[y, x]) > 1:

                    candidateValues = np.zeros([len(vert_coord[y, x])])
                    candidateDepths = np.zeros([len(vert_coord[y, x])])

                    for z in range(0, len(vert_coord[y, x])):
                        # On retrouve l'index de la layer
                        array = np.asarray(indexes_z)
                        index_layer = (np.abs(array - vert_coord[y, x][z])).argmin()

                        if self.is_sigma_coordinate(type="source"):
                            candidateDepths[z] = z_axis[index_layer, y, x]
                        else:
                            candidateDepths[z] = z_axis[index_layer]

                        candidateValues[z] = layers[index_layer, y, x]

                    data[y, x] = vertical_interpolation(candidateDepths, targetDepth, candidateValues,
                                                        LevelCoverage.VERTICAL_INTERPOLATION_METHOD)

        if self.horizontal_resampling:
            data = resample_2d_to_grid(self.read_axis_x(type="source_global", with_overlap=True),
                                       self.read_axis_y(type="source_global", with_overlap=True),
                                       self.read_axis_x(type="target", with_overlap=True),
                                       self.read_axis_y(type="target", with_overlap=True),
                                       data,
                                       Coverage.HORIZONTAL_INTERPOLATION_METHOD)

        return data[self.map_mpi[self.rank]["dst_local_y"], self.map_mpi[self.rank]["dst_local_x"]]


    def read_variable_baroclinic_sea_water_velocity_at_time_and_depth(self,time,depth):
        """Retourne les composantes u,v du courant à la date souhaitée et au niveau souhaité sur toute la couverture horizontale.
    @type time: datetime ou l'index
    @param time: date souhaitée
    @type depth: profondeur en mètre (float) ou index (integer)
    @param depth: profondeur souhaitée. Si le z est un entier, on considère qu'il s'agit de l'index,
    si c'est un flottant on considère qu'il s'agit d'une profondeur
    @return: un tableau en deux dimensions [u_comp,v_comp] contenant chacun deux dimensions [y,x]."""

        index_t = self.find_time_index(time);
        tmp = self.find_level_index(depth);
        vert_coord = tmp[0]
        indexes_z = tmp[1]

        xmax = self.get_x_size(type="source", with_overlap=True)
        ymax = self.get_y_size(type="source", with_overlap=True)

        layers = np.zeros([np.shape(indexes_z)[0], ymax, xmax])
        layers[::] = np.NAN

        data = np.zeros([2,ymax, xmax])
        data[:] = np.NAN
        targetDepth = [depth]

        mask_t = self.reader.read_variable_2D_sea_binary_mask(
            self.map_mpi[self.rank]["src_global_x_overlap"].start,
            self.map_mpi[self.rank]["src_global_x_overlap"].stop,
            self.map_mpi[self.rank]["src_global_y_overlap"].start,
            self.map_mpi[self.rank]["src_global_y_overlap"].stop)

        for z in range(0, len(indexes_z)):
            layers[z] = self.reader.read_variable_baroclinic_sea_water_velocity_at_time_and_depth(
                self.map_mpi[self.rank]["src_global_t"].start + index_t, indexes_z[z],
                self.map_mpi[self.rank]["src_global_x_overlap"].start,
                self.map_mpi[self.rank]["src_global_x_overlap"].stop,
                self.map_mpi[self.rank]["src_global_y_overlap"].start,
                self.map_mpi[self.rank]["src_global_y_overlap"].stop)

        for y in range(0, self.get_y_size(type="source", with_overlap=True)):
            for x in range(0, self.get_x_size(type="source", with_overlap=True)):

                if mask_t[y, x] == 1:

                    if vert_coord[y, x] is not None and len(vert_coord[y, x]) == 1:

                        # Il n'y a qu'une seule couche de sélectionner donc pas d'interpolation possible
                        # On retrouve l'index de la layer
                        array = np.asarray(indexes_z)
                        index_layer = (np.abs(array - vert_coord[y, x][0])).argmin()

                        data[0][y, x] = layers[index_layer,0, y, x]
                        data[1][y, x] = layers[index_layer,1, y, x]

                    elif vert_coord[y, x] is not None and len(vert_coord[y, x]) > 1:

                        candidateValues = np.zeros([2,len(vert_coord[y, x])])
                        candidateDepths = np.zeros([len(vert_coord[y, x])])

                        for z in range(0, len(vert_coord[y, x])):
                            # On retrouve l'index de la layer
                            array = np.asarray(indexes_z)
                            index_layer = (np.abs(array - vert_coord[y, x][z])).argmin()

                            candidateDepths[0][z] = self.source_global_axis_z[index_layer,0, y, x]
                            candidateValues[1][z] = layers[index_layer,1, y, x]

                        data[0][y, x] = vertical_interpolation(candidateDepths, targetDepth, candidateValues[0],
                                                            LevelCoverage.VERTICAL_INTERPOLATION_METHOD)

                        data[1][y, x] = vertical_interpolation(candidateDepths, targetDepth, candidateValues[1],
                                                               LevelCoverage.VERTICAL_INTERPOLATION_METHOD)

        if self.horizontal_resampling:
            data[0] = resample_2d_to_grid(self.read_axis_x(type="source_global", with_overlap=True),
                                          self.read_axis_y(type="source_global", with_overlap=True),
                                          self.read_axis_x(type="target", with_overlap=True),
                                          self.read_axis_y(type="target", with_overlap=True),
                                          data[0],
                                          Coverage.HORIZONTAL_INTERPOLATION_METHOD)

            data[1] = resample_2d_to_grid(self.read_axis_x(type="source_global", with_overlap=True),
                                          self.read_axis_y(type="source_global", with_overlap=True),
                                          self.read_axis_x(type="target", with_overlap=True),
                                          self.read_axis_y(type="target", with_overlap=True),
                                          data[1],
                                          Coverage.HORIZONTAL_INTERPOLATION_METHOD)

        return [data[0][self.map_mpi[self.rank]["dst_local_y"], self.map_mpi[self.rank]["dst_local_x"]],
                data[1][self.map_mpi[self.rank]["dst_local_y"], self.map_mpi[self.rank]["dst_local_x"]]]
