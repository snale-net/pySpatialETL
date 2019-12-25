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
from coverage.Coverage import  Coverage
from coverage.LevelCoverage import LevelCoverage
from coverage.TimeCoverage import TimeCoverage
from coverage.operator.interpolator.InterpolatorCore import resample_2d_to_grid
from coverage.operator.interpolator.InterpolatorCore import vertical_interpolation
import numpy as np

class TimeLevelCoverage(LevelCoverage,TimeCoverage):
    """La classe TimeLevelCoverage est une extension de la classe Coverage, LevelCoverage, TimeCoverage.
Elle rajoute les dimensions temporelle et verticale à la couverture horizontale classique.
    """
    def __init__(self, myReader,bbox=None,resolution_x=None,resolution_y=None,zbox=None,resolution_z=None):

        LevelCoverage.__init__(self,myReader,bbox=bbox,resolution_x=resolution_x,resolution_y=resolution_y,zbox=zbox,resolution_z=resolution_z);
        TimeCoverage.__init__(self,myReader,bbox=bbox,resolution_x=resolution_x,resolution_y=resolution_y);

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

        xmax = self.get_x_size()
        ymax = self.get_y_size()
        layers = np.zeros([np.shape(indexes_z)[0], ymax, xmax])
        layers[::] = np.NAN

        results = np.zeros([ymax, xmax])
        results[:] = np.NAN

        targetDepth = [depth]

        for z in range(0, len(indexes_z)):
            layers[z] = self.reader.read_variable_sea_water_temperature_at_time_and_depth(index_t, indexes_z[z])

        for y in range(0, ymax):
            for x in range(0, xmax):

                if len(vert_coord[y, x]) == 1:
                    # Il n'y a qu'une seule couche de sélectionner donc pas d'interpolation possible

                    # On retrouve l'index de la layer
                    array = np.asarray(indexes_z)
                    index_layer = (np.abs(array - vert_coord[y, x][0])).argmin()

                    results[y, x] = layers[index_layer, y, x]

                elif len(vert_coord[y, x]) > 1:

                    candidateValues = np.zeros([len(vert_coord[y, x])])
                    candidateDepths = np.zeros([len(vert_coord[y, x])])

                    for z in range(0, len(vert_coord[y, x])):
                        # On retrouve l'index de la layer
                        array = np.asarray(indexes_z)
                        index_layer = (np.abs(array - vert_coord[y, x][z])).argmin()

                        candidateDepths[z] = self.levels[index_layer, y, x]
                        candidateValues[z] = layers[index_layer, y, x]

                    results[y, x] = vertical_interpolation(candidateDepths, targetDepth, candidateValues,LevelCoverage.VERTICAL_INTERPOLATION_METHOD)

        return results

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

        xmax = self.get_x_size(type="source",with_overlap=True)
        ymax = self.get_y_size(type="source",with_overlap=True)

        layers = np.zeros([np.shape(indexes_z)[0], ymax, xmax])
        layers[::] = np.NAN

        data = np.zeros([ymax, xmax])
        data[:] = np.NAN

        targetDepth = [depth]

        for z in range(0, len(indexes_z)):
            layers[z] = self.reader.read_variable_sea_water_salinity_at_time_and_depth(
                self.map_mpi[self.rank]["src_global_t"].start+index_t, indexes_z[z],
                self.map_mpi[self.rank]["src_global_x_overlap"].start,
                self.map_mpi[self.rank]["src_global_x_overlap"].stop,
                self.map_mpi[self.rank]["src_global_y_overlap"].start,
                self.map_mpi[self.rank]["src_global_y_overlap"].stop)

        for y in range(0, self.get_y_size(type="source",with_overlap=True)):
            for x in range(0, self.get_x_size(type="source",with_overlap=True)):

                print(vert_coord[y, x] )

                if vert_coord[y, x] is not None and len(vert_coord[y, x]) == 1:
                    # Il n'y a qu'une seule couche de sélectionner donc pas d'interpolation possible

                    # On retrouve l'index de la layer
                    array = np.asarray(indexes_z)
                    index_layer = (np.abs(array - vert_coord[y, x][0])).argmin()

                    data[y, x] = layers[index_layer, y, x]

                elif vert_coord[y, x] is not None and len(vert_coord[y, x]) > 1:

                    candidateValues = np.zeros([len(vert_coord[y, x])])
                    candidateDepths = np.zeros([len(vert_coord[y, x])])

                    for z in range(0, len(vert_coord[y, x])):
                        # On retrouve l'index de la layer
                        array = np.asarray(indexes_z)
                        index_layer = (np.abs(array - vert_coord[y, x][z])).argmin()

                        candidateDepths[z] = self.source_global_axis_z[index_layer, y, x]
                        candidateValues[z] = layers[index_layer, y, x]

                    data[y, x] = vertical_interpolation(candidateDepths, targetDepth, candidateValues,LevelCoverage.VERTICAL_INTERPOLATION_METHOD)

        if self.horizontal_resampling:

            data = resample_2d_to_grid(self.read_axis_x(type="source", with_overlap=True),
                                       self.read_axis_y(type="source", with_overlap=True),
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

        xmax = self.get_x_size()
        ymax = self.get_y_size()
        layers = np.zeros([np.shape(indexes_z)[0],2, ymax, xmax])
        layers[::] = np.NAN

        results = np.zeros([2,ymax, xmax])
        results[:] = np.NAN

        targetDepth = [depth]

        for z in range(0, len(indexes_z)):
            layers[z] = self.reader.read_variable_baroclinic_sea_water_velocity_at_time_and_depth(index_t, indexes_z[z])

        for y in range(0, ymax):
            for x in range(0, xmax):

                if len(vert_coord[y, x]) == 1:
                    # Il n'y a qu'une seule couche de sélectionner donc pas d'interpolation possible

                    # On retrouve l'index de la layer
                    array = np.asarray(indexes_z)
                    index_layer = (np.abs(array - vert_coord[y, x][0])).argmin()

                    results[0][y, x] = layers[index_layer,0, y, x]
                    results[1][y, x] = layers[index_layer,1, y, x]

                elif len(vert_coord[y, x]) > 1:

                    candidateValues = np.zeros([2,len(vert_coord[y, x])])
                    candidateDepths = np.zeros([len(vert_coord[y, x])])

                    for z in range(0, len(vert_coord[y, x])):

                        # On retrouve l'index de la layer
                        array = np.asarray(indexes_z)
                        index_layer = (np.abs(array - vert_coord[y, x][z])).argmin()

                        if self.is_sigma_coordinate():
                            candidateDepths[z] = self.levels[index_layer, y, x]
                        else:
                            candidateDepths[z] = self.levels[index_layer]


                        candidateValues[0][z] = layers[index_layer,0, y, x]
                        candidateValues[1][z] = layers[index_layer,1, y, x]

                    results[0][y, x] = vertical_interpolation(candidateDepths, targetDepth, candidateValues[0],LevelCoverage.VERTICAL_INTERPOLATION_METHOD)
                    results[1][y, x] = vertical_interpolation(candidateDepths, targetDepth,  candidateValues[1],LevelCoverage.VERTICAL_INTERPOLATION_METHOD)

        return results



       

            
        
        
    

