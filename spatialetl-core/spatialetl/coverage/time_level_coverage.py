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
import concurrent
import inspect
import os
from concurrent.futures import ProcessPoolExecutor

import numpy as np

from spatialetl.coverage.coverage import Coverage
from spatialetl.coverage.level_coverage import LevelCoverage
from spatialetl.coverage.time_coverage import TimeCoverage
from spatialetl.operator.interpolator.interpolator_core import resample_2d_to_grid, resample_faster_2d_to_grid
from spatialetl.operator.interpolator.interpolator_core import vertical_interpolation
from spatialetl.utils.logger import logging
from spatialetl.utils.timing import timing


class TimeLevelCoverage(LevelCoverage, TimeCoverage):
    """
    La classe TimeLevelCoverage est une extension de la classe Coverage, LevelCoverage, TimeCoverage.
    Elle rajoute les dimensions temporelle et verticale à la couverture horizontale classique.
    """

    def __init__(self, reader, bbox=None, resolution_x=None, resolution_y=None, zbox=None, resolution_z=None,
                 start_time=None, end_time=None, freq=None, nb_thread:int=os.cpu_count()-1):

        LevelCoverage.__init__(self, reader, bbox=bbox, resolution_x=resolution_x, resolution_y=resolution_y,
                               zbox=zbox, resolution_z=resolution_z, nb_thread=nb_thread);
        TimeCoverage.__init__(self, reader, bbox=bbox, resolution_x=resolution_x, resolution_y=resolution_y,
                              start_time=start_time, end_time=end_time, freq=freq, nb_thread=nb_thread);

        self.data_temp = np.zeros(
            [2, self.get_y_size(type="source", with_overlap=True), self.get_x_size(type="source", with_overlap=True)])
        self.layers_temp = np.zeros(
            [self.get_z_size(type="source"), 2, self.get_y_size(type="source", with_overlap=True),
             self.get_x_size(type="source", with_overlap=True)])

        if self.horizontal_resampling and self.rank == 0:
            logging.info(
                '[horizontal_interpolation] Source grid size : (' + str(self.source_global_x_size) + ", " + str(
                    self.source_global_y_size) + ")")
            logging.info(
                '[horizontal_interpolation] Target grid size : (' + str(self.target_global_x_size) + ", " + str(
                    self.target_global_y_size) + ")")

            self.compute_weight()


    def __read_variable(self, function_name, time, depth):

        fn = getattr(self.reader, function_name)

        index_t = self.find_time_index(time);
        vert_coord, indexes_z = self.find_level_index(depth);
        self.layers_temp[::] = np.nan
        self.data_temp[::] = np.nan
        targetDepth = [depth]

        for z in range(0, len(indexes_z)):
            self.layers_temp[z] = fn(
                self.parallel_map[self.rank]["src_global_t"].start + index_t, indexes_z[z],
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

        if self.horizontal_resampling:

            # We use multithreading to compute resampling in parallel
            local_data = np.zeros([self.get_y_size(), self.get_x_size()])
            local_data[:] = np.nan

            with ProcessPoolExecutor(max_workers=self.threads_number) as executor:
                futures = []
                for current_thread in range(0, executor._max_workers):
                    futures.append(executor.submit(resample_faster_2d_to_grid,
                                                   self.tri[self.rank][current_thread],
                                                   self.read_thread_axis_x(type="target", with_overlap=True,
                                                                           current_thread=current_thread),
                                                   self.read_thread_axis_y(type="target", with_overlap=True,
                                                                           current_thread=current_thread),
                                                   self.data_temp[0],
                                                   Coverage.HORIZONTAL_INTERPOLATION_METHOD,
                                                   current_thread))

                futures, _ = concurrent.futures.wait(futures)

                for f in futures:
                    current_thread, data = f.result()
                    local_data[self.parallel_map[self.rank]["threads_map"][current_thread]["dst_global_y"],
                    self.parallel_map[self.rank]["threads_map"][current_thread]["dst_global_x"]] = data[
                        self.parallel_map[self.rank]["threads_map"][current_thread]["dst_local_y"],
                        self.parallel_map[self.rank]["threads_map"][current_thread]["dst_local_x"]]

            return local_data

        else:
            return self.data_temp[0, self.parallel_map[self.rank]["dst_local_y"], self.parallel_map[self.rank]["dst_local_x"]]


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
        return self.__read_variable(inspect.stack()[0][3],time=time, depth=depth)




    @timing
    def read_variable_sea_water_salinity_at_time_and_depth(self, time, depth):
        """Retourne la salinité à la date souhaitée et au niveau souhaité sur toute la couverture horizontale.
    @type time: datetime ou l'index
    @param time: date souhaitée
    @type depth: profondeur en mètre (float) ou index (integer)
    @param depth: profondeur souhaitée. Si le z est un entier, on considère qu'il s'agit de l'index,
    si c'est un flottant on considère qu'il s'agit d'une profondeur
    @return: un tableau en deux dimensions [y,x]."""

        return self.__read_variable(inspect.stack()[0][3],time=time, depth=depth)

    def read_variable_baroclinic_sea_water_velocity_at_time_and_depth(self, time, depth):
        """Retourne les composantes u,v du courant à la date souhaitée et au niveau souhaité sur toute la couverture horizontale.
    @type time: datetime ou l'index
    @param time: date souhaitée
    @type depth: profondeur en mètre (float) ou index (integer)
    @param depth: profondeur souhaitée. Si le z est un entier, on considère qu'il s'agit de l'index,
    si c'est un flottant on considère qu'il s'agit d'une profondeur
    @return: un tableau en deux dimensions [u_comp,v_comp] contenant chacun deux dimensions [y,x]."""

        index_t = self.find_time_index(time);
        vert_coord, indexes_z = self.find_level_index(depth);
        self.layers_temp[::] = np.nan
        self.data_temp[:] = np.nan
        targetDepth = [depth]

        for z in range(0, len(indexes_z)):
            self.layers_temp[z] = self.reader.read_variable_baroclinic_sea_water_velocity_at_time_and_depth(
                self.parallel_map[self.rank]["src_global_t"].start + index_t, indexes_z[z],
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
                self.data_temp[1, y, x] = self.layers_temp[index_layer, 1, y, x]
            else:
                candidateValues = np.zeros([2, len(vert_coord[y, x])])
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

                    candidateValues[0, z] = self.layers_temp[index_layer, 0, y, x]
                    candidateValues[1, z] = self.layers_temp[index_layer, 1, y, x]

                self.data_temp[0, y, x] = vertical_interpolation(candidateDepths, targetDepth, candidateValues[0],
                                                                 LevelCoverage.VERTICAL_INTERPOLATION_METHOD)

                self.data_temp[1, y, x] = vertical_interpolation(candidateDepths, targetDepth, candidateValues[1],
                                                                 LevelCoverage.VERTICAL_INTERPOLATION_METHOD)

        if self.horizontal_resampling:
            return resample_2d_to_grid(self.read_axis_x(type="source", with_overlap=True),
                                       self.read_axis_y(type="source", with_overlap=True),
                                       self.read_axis_x(type="target", with_overlap=True),
                                       self.read_axis_y(type="target", with_overlap=True),
                                       self.data_temp[0],
                                       Coverage.HORIZONTAL_INTERPOLATION_METHOD)[
                self.parallel_map[self.rank]["dst_local_y"], self.parallel_map[self.rank]["dst_local_x"]], \
            resample_2d_to_grid(self.read_axis_x(type="source", with_overlap=True),
                                self.read_axis_y(type="source", with_overlap=True),
                                self.read_axis_x(type="target", with_overlap=True),
                                self.read_axis_y(type="target", with_overlap=True),
                                self.data_temp[1],
                                Coverage.HORIZONTAL_INTERPOLATION_METHOD)[
                self.parallel_map[self.rank]["dst_local_y"], self.parallel_map[self.rank]["dst_local_x"]]

        return self.data_temp[0, self.parallel_map[self.rank]["dst_local_y"], self.parallel_map[self.rank]["dst_local_x"]], \
        self.data_temp[1, self.parallel_map[self.rank]["dst_local_y"], self.parallel_map[self.rank]["dst_local_x"]]
