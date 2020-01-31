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
from datetime import timedelta
from datetime import datetime
import cftime
from numpy import int,int32,int64
import math
import numpy as np
from array_split import shape_split
from spatialetl.operator.interpolator.InterpolatorCore import resample_2d_to_grid
import logging
import pandas

class TimeCoverage(Coverage):
    """La classe TimeCoverage est une extension de la classe Coverage.
Elle rajoute une dimension temporelle à la couverture horizontale classique.
    """
    
    TIME_DATUM = datetime(1970, 1, 1)    
    TIME_DELTA = timedelta(minutes = 15)
    TIME_OVERLAPING_SIZE = 0

    def __init__(self, myReader,bbox=None,resolution_x=None,resolution_y=None,start_time=None,end_time=None,freq=None):

        Coverage.__init__(self, myReader, bbox=bbox, resolution_x=resolution_x, resolution_y=resolution_y);

        self.source_global_t_size = self.reader.get_t_size()
        self.source_global_axis_t = self.reader.read_axis_t(0,self.source_global_t_size,0);

        self.temporal_resampling = False
        tmin = 0
        tmax = self.source_global_t_size
        zero_delta = timedelta(minutes=00)

        if start_time is not None:

            if type(start_time) == datetime or type(start_time) == cftime._cftime.real_datetime:
                time = start_time
            elif type(start_time) == str:
                try:
                    time = datetime.strptime(start_time, '%Y-%m-%d %H:%M:%S')
                except ValueError as ex:
                    raise ValueError("start_time is not well formated " + str(ex))
            else:
                raise ValueError("start_time have to be string or datetime. Found " + str(type(start_time)))

            nearest_t_index = (np.abs(self.source_global_axis_t - time)).argmin()

            if time - self.source_global_axis_t[nearest_t_index] == zero_delta or abs(time - self.source_global_axis_t[nearest_t_index]) < TimeCoverage.TIME_DELTA:
                tmin = nearest_t_index
            else:
                raise ValueError(str(time) + " not found. Maybe the TimeCoverage.TIME_DELTA (" + str(
                TimeCoverage.TIME_DELTA) + ") is too small or the date is out the range.")

        if end_time is not None:

            self.temporal_resampling = True

            if type(end_time) == datetime or type(end_time) == cftime._cftime.real_datetime:
                time = end_time
            elif type(end_time) == str:
                try:
                    time = datetime.strptime(end_time, '%Y-%m-%d %H:%M:%S')
                except ValueError as ex:
                    raise ValueError("end_time is not well formated "+str(ex))
            else:
                raise ValueError("end_time have to be string or datetime. Found " + str(type(end_time)))

            nearest_t_index = (np.abs(self.source_global_axis_t - time)).argmin()

            if time - self.source_global_axis_t[nearest_t_index] == zero_delta or abs(time - self.source_global_axis_t[nearest_t_index]) < TimeCoverage.TIME_DELTA:
                tmax = nearest_t_index +1
            else:
                raise ValueError(str(time) + " not found. Maybe the TimeCoverage.TIME_DELTA (" + str(
                    TimeCoverage.TIME_DELTA) + ") is too small or the date is out the range.")

        if freq is not None:
            self.target_global_axis_t = pandas.date_range(start=self.source_global_axis_t[tmin],
                                                  end=self.source_global_axis_t[tmax-1], freq=freq).to_pydatetime();
            self.target_global_t_size = np.shape(self.target_global_axis_t)[0]
        else:
            self.target_global_axis_t = self.source_global_axis_t[tmin:tmax]
            self.target_global_t_size = tmax - tmin

        self.create_mpi_map()
        self.update_mpi_map()

        if type(self) == TimeCoverage and self.horizontal_resampling and self.rank == 0:
                logging.info(
                    '[horizontal_interpolation] Source grid size : (' + str(self.source_global_x_size) + ", " + str(
                        self.source_global_y_size) + ")")
                logging.info(
                    '[horizontal_interpolation] Target grid size : (' + str(self.target_global_x_size) + ", " + str(
                        self.target_global_y_size) + ")")

        #self.map_mpi[self.rank]["src_global_t"] = np.s_[tmin:tmax]
        #self.map_mpi[self.rank]["src_global_overlap"] = np.s_[tmin:tmax]
        #self.map_mpi[self.rank]["src_local_t"] = np.s_[0:self.source_global_t_size]

        if self.rank == 0:
            logging.debug("MPI map:")
        for key in self.map_mpi[self.rank]:
            logging.debug("Proc n°" + str(self.rank) + " " + str(key) + "=" + str(self.map_mpi[self.rank][key]))
        logging.debug("---------")

    def create_mpi_map(self):

        self.map_mpi = np.empty(self.size, dtype=object)
        target_sample = (self.target_global_t_size, self.target_global_y_size, self.target_global_x_size)

        # Découpage des axes
        # if self.horizontal_resampling:
        #     # Découpage sur le temps uniquemenent
        #     target_slices = shape_split(target_sample, self.size, axis=[0, 1, 1])
        #     # Si on feet pas le nombre de proc
        #     if len( target_slices.flatten()) != self.size:
        #         target_slices = shape_split(target_sample, self.size, axis=[0, 0, 0])
        # else:

        target_slices = shape_split(target_sample, self.size, axis=[0, 0, 0])

        slice_index = 0
        for slyce in target_slices.flatten():
            slice = tuple(slyce)

            map = {}
            # Grille source
            map["dst_global_t"] = slice[0]
            map["dst_global_x"] = slice[2]
            map["dst_global_y"] = slice[1]

            map["dst_local_t_size"] = map["dst_global_t"].stop - map["dst_global_t"].start
            map["dst_local_x_size"] = map["dst_global_x"].stop - map["dst_global_x"].start
            map["dst_local_y_size"] = map["dst_global_y"].stop - map["dst_global_y"].start

            dst_global_t_min_overlap = max(0, map["dst_global_t"].start - TimeCoverage.TIME_OVERLAPING_SIZE)
            dst_global_t_max_overlap = min(self.target_global_t_size, map["dst_global_t"].stop + TimeCoverage.TIME_OVERLAPING_SIZE)
            map["dst_global_t_overlap"] = np.s_[dst_global_t_min_overlap:dst_global_t_max_overlap]

            dst_global_x_min_overlap = max(0,map["dst_global_x"].start - Coverage.HORIZONTAL_OVERLAPING_SIZE)
            dst_global_x_max_overlap = min(self.target_global_x_size, map["dst_global_x"].stop + Coverage.HORIZONTAL_OVERLAPING_SIZE)
            map["dst_global_x_overlap"] = np.s_[dst_global_x_min_overlap:dst_global_x_max_overlap]

            dst_global_y_min_overlap = max(0,map["dst_global_y"].start - Coverage.HORIZONTAL_OVERLAPING_SIZE)
            dst_global_y_max_overlap = min(self.target_global_y_size, map["dst_global_y"].stop + Coverage.HORIZONTAL_OVERLAPING_SIZE)
            map["dst_global_y_overlap"] = np.s_[dst_global_y_min_overlap:dst_global_y_max_overlap]

            map["dst_global_t_size_overlap"] = map["dst_global_t_overlap"].stop - map["dst_global_t_overlap"].start
            map["dst_global_x_size_overlap"] = map["dst_global_x_overlap"].stop - map["dst_global_x_overlap"].start
            map["dst_global_y_size_overlap"] = map["dst_global_y_overlap"].stop - map["dst_global_y_overlap"].start

            dst_t_min = TimeCoverage.TIME_OVERLAPING_SIZE
            dst_t_max = map["dst_global_t_size_overlap"] - TimeCoverage.TIME_OVERLAPING_SIZE
            dst_x_min = Coverage.HORIZONTAL_OVERLAPING_SIZE
            dst_x_max = map["dst_global_x_size_overlap"] - Coverage.HORIZONTAL_OVERLAPING_SIZE
            dst_y_min = Coverage.HORIZONTAL_OVERLAPING_SIZE
            dst_y_max = map["dst_global_y_size_overlap"] - Coverage.HORIZONTAL_OVERLAPING_SIZE

            if map["dst_global_t"].start == 0:
                dst_t_min = 0

            if map["dst_global_t"].stop == self.target_global_t_size:
                dst_t_max = map["dst_global_t_size_overlap"]

            if map["dst_global_x"].start == 0:
                dst_x_min = 0

            if map["dst_global_x"].stop == self.target_global_x_size:
                dst_x_max = map["dst_global_x_size_overlap"]

            if map["dst_global_y"].start == 0:
                dst_y_min = 0

            if map["dst_global_y"].stop == self.target_global_y_size:
                dst_y_max = map["dst_global_y_size_overlap"]

            map["dst_local_t"] = np.s_[dst_t_min:dst_t_max]
            map["dst_local_x"] = np.s_[dst_x_min:dst_x_max]
            map["dst_local_y"] = np.s_[dst_y_min:dst_y_max]

            # Source grille
            map["src_global_t"] = map["dst_global_t"]
            map["src_global_x"] = map["dst_global_x"]
            map["src_global_y"] = map["dst_global_y"]

            map["src_global_t_overlap"] = map["dst_global_t_overlap"]
            map["src_global_x_overlap"] = map["dst_global_x_overlap"]
            map["src_global_y_overlap"] = map["dst_global_y_overlap"]

            map["src_local_t"] = map["dst_local_t"]
            map["src_local_x"] = map["dst_local_x"]
            map["src_local_y"] = map["dst_local_y"]

            map["src_local_t_size"] = map["dst_local_t_size"]
            map["src_local_x_size"] = map["dst_local_x_size"]
            map["src_local_y_size"] = map["dst_local_y_size"]

            map["src_local_t_size_overlap"] = map["dst_global_t_size_overlap"]
            map["src_local_x_size_overlap"] = map["dst_global_x_size_overlap"]
            map["src_local_y_size_overlap"] = map["dst_global_y_size_overlap"]

            self.map_mpi[slice_index] = map

            slice_index = slice_index + 1
   
    # Axis
    def find_time_index(self,t,method="quick"):
        """Retourne l'index de la date la plus proche à TIME_DELTA_MIN prêt.
    @type t: datetime ou int
    @param t: date souhaitée ou l'index de la date souhaitée
    @return:  l'index de la date la plus proche à TIME_DELTA_MIN prêt ou une erreur si aucune date n'a pu être trouvée."""

        if type(t) == int or type(t) == int32 or type(t) == int64:

            if t < 0 or t >= self.get_t_size():
                raise ValueError("Time index have to range between 0 and " + str(
                    self.get_t_size() - 1) + ". Actually Time index = " + str(t))

            return t;

        if type(t) == datetime or type(t) == cftime._cftime.real_datetime:

            zero_delta = timedelta(minutes=00)
            array = np.asarray(self.target_global_axis_t[self.map_mpi[self.rank]["src_global_t"]])

            #logging.debug("[TimeCoverage][find_time_index()] Looking for : " + str(t))
            #logging.debug("[TimeCovergae][find_time_index()] Proc [" + str(self.rank) + "] - datetime candidates are : " + str(array))

            if method == "classic":

                for i in range(0,self.get_t_size()):
                    if t - array[i] == zero_delta or abs(t - array[i]) < TimeCoverage.TIME_DELTA:
                        return i

            elif method == "quick":

                nearest_t_index = (np.abs(array - t)).argmin()

                if t - array[nearest_t_index] == zero_delta or abs(t - array[nearest_t_index]) < TimeCoverage.TIME_DELTA:
                    return nearest_t_index
            else:
                raise NotImplementedError("Method " + str(method) + " is not implemented for regular grid.")

            raise ValueError("Proc n°"+str(self.rank)+" doesn't find '"+str(t)+"'. Maybe the TimeCoverage.TIME_DELTA ("+ str(TimeCoverage.TIME_DELTA)+") is too small or the date is out the range.")

        else:
            raise ValueError(""+str(t)+" have to be an integer or a datetime. Current type: "+str(type(t)))
    
    def read_axis_t(self,type="target",with_overlap=False,timestamp=0):
        """Retourne les valeurs de l'axe t.
    @param timestamp: égale 1 si le temps est souhaité en timestamp depuis TIME_DATUM.
    @return:  un tableau à une dimensions [z] au format datetime ou timestamp si timestamp=1."""
        if type == "target_global":
            return self.target_global_axis_t

        elif type == "source_global":
            return self.reader.read_axis_t(self.map_mpi[self.rank]["src_global_t"].start,
                                           self.map_mpi[self.rank]["src_global_t"].stop,timestamp)

        elif type == "source":
            return self.reader.read_axis_t(self.map_mpi[self.rank]["src_local_t"].start,
                                           self.map_mpi[self.rank]["src_local_t"].stop,timestamp)

        elif type == "target" and with_overlap is True:
                return self.target_global_axis_t[self.map_mpi[self.rank]["dst_global_t_overlap"]]
        else:
                return self.target_global_axis_t[self.map_mpi[self.rank]["dst_global_t"]]

    def get_t_size(self,type="target",with_overlap=False):
        if type == "target_global":
            return self.target_global_t_size
        elif type == "source_global":
            return self.source_global_t_size
        elif type == "source":
            return self.map_mpi[self.rank]["src_local_t_size"]
        elif type == "target" and with_overlap is True:
            return self.map_mpi[self.rank]["dst_local_t_size_overlap"]
        else:
            return self.map_mpi[self.rank]["dst_local_t_size"]
    
    # Variables
    def read_variable_2D_sea_binary_mask_at_time(self, t):
        """Retourne le masque à la date souhaitée sur toute la couverture horizontale.
    @type t: datetime ou l'index
    @param t: date souhaitée
    @return: un tableau en deux dimensions [y,x]."""

        index_t = self.find_time_index(t);

        data = self.reader.read_variable_2D_sea_binary_mask_at_time(
            self.map_mpi[self.rank]["src_global_t"].start + index_t,
            self.map_mpi[self.rank]["src_global_x_overlap"].start,
            self.map_mpi[self.rank]["src_global_x_overlap"].stop,
            self.map_mpi[self.rank]["src_global_y_overlap"].start,
            self.map_mpi[self.rank]["src_global_y_overlap"].stop)

        if self.horizontal_resampling:
            data = resample_2d_to_grid(self.read_axis_x(type="source_global", with_overlap=True),
                                       self.read_axis_y(type="source_global", with_overlap=True),
                                       self.read_axis_x(type="target", with_overlap=True),
                                       self.read_axis_y(type="target", with_overlap=True),
                                       data,
                                       Coverage.HORIZONTAL_INTERPOLATION_METHOD)

        return data[self.map_mpi[self.rank]["dst_local_y"], self.map_mpi[self.rank]["dst_local_x"]]

    def read_variable_2D_wet_binary_mask_at_time(self, t):
        """Retourne le masque à la date souhaitée sur toute la couverture horizontale.
    @type t: datetime ou l'index
    @param t: date souhaitée
    @return: un tableau en deux dimensions [y,x]."""

        index_t = self.find_time_index(t);

        data = self.reader.read_variable_2D_wet_binary_mask_at_time(
            self.map_mpi[self.rank]["src_global_t"].start + index_t,
            self.map_mpi[self.rank]["src_global_x_overlap"].start,
            self.map_mpi[self.rank]["src_global_x_overlap"].stop,
            self.map_mpi[self.rank]["src_global_y_overlap"].start,
            self.map_mpi[self.rank]["src_global_y_overlap"].stop)

        if self.horizontal_resampling:
            data = resample_2d_to_grid(self.read_axis_x(type="source_global", with_overlap=True),
                                       self.read_axis_y(type="source_global", with_overlap=True),
                                       self.read_axis_x(type="target", with_overlap=True),
                                       self.read_axis_y(type="target", with_overlap=True),
                                       data,
                                       Coverage.HORIZONTAL_INTERPOLATION_METHOD)

        return data[self.map_mpi[self.rank]["dst_local_y"], self.map_mpi[self.rank]["dst_local_x"]]

    #################
    # HYDRO
    # Sea Surface
    #################
    def read_variable_sea_surface_height_above_mean_sea_level_at_time(self,t):

        index_t = self.find_time_index(t);

        data = self.reader.read_variable_sea_surface_height_above_mean_sea_level_at_time(
            self.map_mpi[self.rank]["src_global_t"].start+index_t,
            self.map_mpi[self.rank]["src_global_x_overlap"].start,
            self.map_mpi[self.rank]["src_global_x_overlap"].stop,
            self.map_mpi[self.rank]["src_global_y_overlap"].start,
            self.map_mpi[self.rank]["src_global_y_overlap"].stop)

        if self.horizontal_resampling:

            data = resample_2d_to_grid(self.read_axis_x(type="source_global", with_overlap=True),
                                       self.read_axis_y(type="source_global", with_overlap=True),
                                       self.read_axis_x(type="target", with_overlap=True),
                                       self.read_axis_y(type="target", with_overlap=True),
                                       data,
                                       Coverage.HORIZONTAL_INTERPOLATION_METHOD)

        return data[self.map_mpi[self.rank]["dst_local_y"],self.map_mpi[self.rank]["dst_local_x"]]

    def read_variable_sea_surface_height_above_geoid_at_time(self, t):

        index_t = self.find_time_index(t);

        data = self.reader.read_variable_sea_surface_height_above_geoid_at_time(
            self.map_mpi[self.rank]["src_global_t"].start + index_t,
            self.map_mpi[self.rank]["src_global_x_overlap"].start,
            self.map_mpi[self.rank]["src_global_x_overlap"].stop,
            self.map_mpi[self.rank]["src_global_y_overlap"].start,
            self.map_mpi[self.rank]["src_global_y_overlap"].stop)

        if self.horizontal_resampling:
            data = resample_2d_to_grid(self.read_axis_x(type="source_global", with_overlap=True),
                                       self.read_axis_y(type="source_global", with_overlap=True),
                                       self.read_axis_x(type="target", with_overlap=True),
                                       self.read_axis_y(type="target", with_overlap=True),
                                       data,
                                       Coverage.HORIZONTAL_INTERPOLATION_METHOD)

        return data[self.map_mpi[self.rank]["dst_local_y"], self.map_mpi[self.rank]["dst_local_x"]]

    def read_variable_sea_surface_temperature_at_time(self, t):
        """Retourne la temperature de surface à la date souhaitée
    @type t: datetime ou l'index
    @param t: date souhaitée
    @return: un tableau en deux dimensions [y,x]."""

        index_t = self.find_time_index(t);

        data = self.reader.read_variable_sea_surface_temperature_at_time(
            self.map_mpi[self.rank]["src_global_t"].start + index_t,
            self.map_mpi[self.rank]["src_global_x_overlap"].start,
            self.map_mpi[self.rank]["src_global_x_overlap"].stop,
            self.map_mpi[self.rank]["src_global_y_overlap"].start,
            self.map_mpi[self.rank]["src_global_y_overlap"].stop)

        if self.horizontal_resampling:
            data = resample_2d_to_grid(self.read_axis_x(type="source_global", with_overlap=True),
                                       self.read_axis_y(type="source_global", with_overlap=True),
                                       self.read_axis_x(type="target", with_overlap=True),
                                       self.read_axis_y(type="target", with_overlap=True),
                                       data,
                                       Coverage.HORIZONTAL_INTERPOLATION_METHOD)

        return data[self.map_mpi[self.rank]["dst_local_y"], self.map_mpi[self.rank]["dst_local_x"]]

    def read_variable_sea_surface_salinity_at_time(self, t):
        """Retourne la salinité de surface à la date souhaitée
    @type t: datetime ou l'index
    @param t: date souhaitée
    @return: un tableau en deux dimensions [y,x]."""

        index_t = self.find_time_index(t);

        data = self.reader.read_variable_sea_surface_salinity_at_time(
            self.map_mpi[self.rank]["src_global_t"].start + index_t,
            self.map_mpi[self.rank]["src_global_x_overlap"].start,
            self.map_mpi[self.rank]["src_global_x_overlap"].stop,
            self.map_mpi[self.rank]["src_global_y_overlap"].start,
            self.map_mpi[self.rank]["src_global_y_overlap"].stop)

        if self.horizontal_resampling:
            data = resample_2d_to_grid(self.read_axis_x(type="source_global", with_overlap=True),
                                       self.read_axis_y(type="source_global", with_overlap=True),
                                       self.read_axis_x(type="target", with_overlap=True),
                                       self.read_axis_y(type="target", with_overlap=True),
                                       data,
                                       Coverage.HORIZONTAL_INTERPOLATION_METHOD)

        return data[self.map_mpi[self.rank]["dst_local_y"], self.map_mpi[self.rank]["dst_local_x"]]

    def read_variable_sea_surface_pressure_at_time(self, t):
        """Retourne la pression à la surface de la mer (sea surface pressure) à la date souhaitée sur toute la couverture horizontale.
    @type t: datetime ou l'index
    @param t: date souhaitée
    @return: un tableau en deux dimensions [y,x]."""

        index_t = self.find_time_index(t);

        data = self.reader.read_variable_sea_surface_pressure_at_time(
            self.map_mpi[self.rank]["src_global_t"].start + index_t,
            self.map_mpi[self.rank]["src_global_x_overlap"].start,
            self.map_mpi[self.rank]["src_global_x_overlap"].stop,
            self.map_mpi[self.rank]["src_global_y_overlap"].start,
            self.map_mpi[self.rank]["src_global_y_overlap"].stop)

        if self.horizontal_resampling:
            data = resample_2d_to_grid(self.read_axis_x(type="source_global", with_overlap=True),
                                       self.read_axis_y(type="source_global", with_overlap=True),
                                       self.read_axis_x(type="target", with_overlap=True),
                                       self.read_axis_y(type="target", with_overlap=True),
                                       data,
                                       Coverage.HORIZONTAL_INTERPOLATION_METHOD)

        return data[self.map_mpi[self.rank]["dst_local_y"], self.map_mpi[self.rank]["dst_local_x"]]

    def read_variable_sea_surface_density_at_time(self, t):
        """Retourne la densité de l'eau de surface à la date souhaitée
    @type t: datetime ou l'index
    @param t: date souhaitée
    @return: un tableau en deux dimensions [y,x]."""

        index_t = self.find_time_index(t);

        data = self.reader.read_variable_sea_surface_density_at_time(
            self.map_mpi[self.rank]["src_global_t"].start + index_t,
            self.map_mpi[self.rank]["src_global_x_overlap"].start,
            self.map_mpi[self.rank]["src_global_x_overlap"].stop,
            self.map_mpi[self.rank]["src_global_y_overlap"].start,
            self.map_mpi[self.rank]["src_global_y_overlap"].stop)

        if self.horizontal_resampling:
            data = resample_2d_to_grid(self.read_axis_x(type="source_global", with_overlap=True),
                                       self.read_axis_y(type="source_global", with_overlap=True),
                                       self.read_axis_x(type="target", with_overlap=True),
                                       self.read_axis_y(type="target", with_overlap=True),
                                       data,
                                       Coverage.HORIZONTAL_INTERPOLATION_METHOD)

        return data[self.map_mpi[self.rank]["dst_local_y"], self.map_mpi[self.rank]["dst_local_x"]]

    def read_variable_sea_water_turbidity_at_time(self, t):
        """Retourne la turbidité de l'eau de surface à la date souhaitée
    @type t: datetime ou l'index
    @param t: date souhaitée
    @return: un tableau en deux dimensions [y,x]."""

        index_t = self.find_time_index(t);

        data = self.reader.read_variable_sea_water_turbidity_at_time(
            self.map_mpi[self.rank]["src_global_t"].start + index_t,
            self.map_mpi[self.rank]["src_global_x_overlap"].start,
            self.map_mpi[self.rank]["src_global_x_overlap"].stop,
            self.map_mpi[self.rank]["src_global_y_overlap"].start,
            self.map_mpi[self.rank]["src_global_y_overlap"].stop)

        if self.horizontal_resampling:
            data = resample_2d_to_grid(self.read_axis_x(type="source_global", with_overlap=True),
                                       self.read_axis_y(type="source_global", with_overlap=True),
                                       self.read_axis_x(type="target", with_overlap=True),
                                       self.read_axis_y(type="target", with_overlap=True),
                                       data,
                                       Coverage.HORIZONTAL_INTERPOLATION_METHOD)

        return data[self.map_mpi[self.rank]["dst_local_y"], self.map_mpi[self.rank]["dst_local_x"]]


    def read_variable_sea_water_velocity_at_sea_water_surface_at_time(self, t):
        """Retourne les composantes u,v du courant à la date souhaitée
    @type t: datetime ou l'index
    @param t: date souhaitée
    @return: un tableau en deux dimensions [u_comp,v_comp] contenant chacun deux dimensions [y,x]."""

        index_t = self.find_time_index(t);

        data = self.reader.read_variable_sea_water_velocity_at_sea_water_surface_at_time(
            self.map_mpi[self.rank]["src_global_t"].start + index_t,
            self.map_mpi[self.rank]["src_global_x_overlap"].start,
            self.map_mpi[self.rank]["src_global_x_overlap"].stop,
            self.map_mpi[self.rank]["src_global_y_overlap"].start,
            self.map_mpi[self.rank]["src_global_y_overlap"].stop)

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

    #################
    # HYDRO
    # Ground level
    #################

    def read_variable_sea_water_temperature_at_ground_level_at_time(self, t):
        """Retourne la temperature de surface à la date souhaitée
    @type t: datetime ou l'index
    @param t: date souhaitée
    @return: un tableau en deux dimensions [y,x]."""

        index_t = self.find_time_index(t);

        data = self.reader.read_variable_sea_water_temperature_at_ground_level_at_time(
            self.map_mpi[self.rank]["src_global_t"].start + index_t,
            self.map_mpi[self.rank]["src_global_x_overlap"].start,
            self.map_mpi[self.rank]["src_global_x_overlap"].stop,
            self.map_mpi[self.rank]["src_global_y_overlap"].start,
            self.map_mpi[self.rank]["src_global_y_overlap"].stop)

        if self.horizontal_resampling:
            data = resample_2d_to_grid(self.read_axis_x(type="source_global", with_overlap=True),
                                       self.read_axis_y(type="source_global", with_overlap=True),
                                       self.read_axis_x(type="target", with_overlap=True),
                                       self.read_axis_y(type="target", with_overlap=True),
                                       data,
                                       Coverage.HORIZONTAL_INTERPOLATION_METHOD)

        return data[self.map_mpi[self.rank]["dst_local_y"], self.map_mpi[self.rank]["dst_local_x"]]

    def read_variable_sea_water_salinity_at_ground_level_at_time(self, t):
        """Retourne la salinité de surface à la date souhaitée
    @type t: datetime ou l'index
    @param t: date souhaitée
    @return: un tableau en deux dimensions [y,x]."""

        index_t = self.find_time_index(t);

        data = self.reader.read_variable_sea_water_salinity_at_ground_level_at_time(
            self.map_mpi[self.rank]["src_global_t"].start + index_t,
            self.map_mpi[self.rank]["src_global_x_overlap"].start,
            self.map_mpi[self.rank]["src_global_x_overlap"].stop,
            self.map_mpi[self.rank]["src_global_y_overlap"].start,
            self.map_mpi[self.rank]["src_global_y_overlap"].stop)

        if self.horizontal_resampling:
            data = resample_2d_to_grid(self.read_axis_x(type="source_global", with_overlap=True),
                                       self.read_axis_y(type="source_global", with_overlap=True),
                                       self.read_axis_x(type="target", with_overlap=True),
                                       self.read_axis_y(type="target", with_overlap=True),
                                       data,
                                       Coverage.HORIZONTAL_INTERPOLATION_METHOD)

        return data[self.map_mpi[self.rank]["dst_local_y"], self.map_mpi[self.rank]["dst_local_x"]]

    def read_variable_sea_water_velocity_at_ground_level_at_time(self, t):
        """Retourne les composantes u,v du courant à la date souhaitée
           @type t: datetime ou l'index
           @param t: date souhaitée
           @return: un tableau en deux dimensions [u_comp,v_comp] contenant chacun deux dimensions [y,x]."""

        index_t = self.find_time_index(t);

        data = self.reader.read_variable_sea_water_velocity_at_ground_level_at_time(
            self.map_mpi[self.rank]["src_global_t"].start + index_t,
            self.map_mpi[self.rank]["src_global_x_overlap"].start,
            self.map_mpi[self.rank]["src_global_x_overlap"].stop,
            self.map_mpi[self.rank]["src_global_y_overlap"].start,
            self.map_mpi[self.rank]["src_global_y_overlap"].stop)

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

    #################
    # HYDRO
    # 2D
    #################
    def read_variable_barotropic_sea_water_velocity_at_time(self,t):
        """Retourne les composantes u,v du courant à la date souhaitée
    @type t: datetime ou l'index
    @param t: date souhaitée
    @return: un tableau en deux dimensions [u_comp,v_comp] contenant chacun deux dimensions [y,x]."""

        index_t = self.find_time_index(t);

        data = self.reader.read_variable_barotropic_sea_water_velocity_at_time(
            self.map_mpi[self.rank]["src_global_t"].start + index_t,
            self.map_mpi[self.rank]["src_global_x_overlap"].start,
            self.map_mpi[self.rank]["src_global_x_overlap"].stop,
            self.map_mpi[self.rank]["src_global_y_overlap"].start,
            self.map_mpi[self.rank]["src_global_y_overlap"].stop)

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

        return [data[0][self.map_mpi[self.rank]["dst_local_y"],self.map_mpi[self.rank]["dst_local_x"]],data[1][self.map_mpi[self.rank]["dst_local_y"],self.map_mpi[self.rank]["dst_local_x"]]]

    #################
    # WAVES
    # Sea Surface
    #################
    def read_variable_sea_surface_wave_significant_height_at_time(self,t):
        """Retourne la hauteur significative des vagues à la date souhaitée sur toute la couverture horizontale.
    @type t: datetime ou l'index
    @param t: date souhaitée
    @return: un tableau en deux dimensions [y,x]."""

        index_t = self.find_time_index(t);

        data = self.reader.read_variable_sea_surface_wave_significant_height_at_time(
            self.map_mpi[self.rank]["src_global_t"].start + index_t,
            self.map_mpi[self.rank]["src_global_x_overlap"].start,
            self.map_mpi[self.rank]["src_global_x_overlap"].stop,
            self.map_mpi[self.rank]["src_global_y_overlap"].start,
            self.map_mpi[self.rank]["src_global_y_overlap"].stop)

        if self.horizontal_resampling:
            data = resample_2d_to_grid(self.read_axis_x(type="source_global", with_overlap=True),
                                       self.read_axis_y(type="source_global", with_overlap=True),
                                       self.read_axis_x(type="target", with_overlap=True),
                                       self.read_axis_y(type="target", with_overlap=True),
                                       data,
                                       Coverage.HORIZONTAL_INTERPOLATION_METHOD)

        return data[self.map_mpi[self.rank]["dst_local_y"], self.map_mpi[self.rank]["dst_local_x"]]

    def read_variable_sea_surface_wave_breaking_height_at_time(self, t):
        """Retourne la hauteur de déferlement des vagues à la date souhaitée sur toute la couverture horizontale.
    @type t: datetime ou l'index
    @param t: date souhaitée
    @return: un tableau en deux dimensions [y,x]."""

        index_t = self.find_time_index(t);

        data = self.reader.read_variable_sea_surface_wave_breaking_height_at_time(
            self.map_mpi[self.rank]["src_global_t"].start + index_t,
            self.map_mpi[self.rank]["src_global_x_overlap"].start,
            self.map_mpi[self.rank]["src_global_x_overlap"].stop,
            self.map_mpi[self.rank]["src_global_y_overlap"].start,
            self.map_mpi[self.rank]["src_global_y_overlap"].stop)

        if self.horizontal_resampling:
            data = resample_2d_to_grid(self.read_axis_x(type="source_global", with_overlap=True),
                                       self.read_axis_y(type="source_global", with_overlap=True),
                                       self.read_axis_x(type="target", with_overlap=True),
                                       self.read_axis_y(type="target", with_overlap=True),
                                       data,
                                       Coverage.HORIZONTAL_INTERPOLATION_METHOD)

        return data[self.map_mpi[self.rank]["dst_local_y"], self.map_mpi[self.rank]["dst_local_x"]]
    
    def read_variable_sea_surface_wave_mean_period_at_time(self,t):

        index_t = self.find_time_index(t);

        data = self.reader.read_variable_sea_surface_wave_mean_period_at_time(
            self.map_mpi[self.rank]["src_global_t"].start + index_t,
            self.map_mpi[self.rank]["src_global_x_overlap"].start,
            self.map_mpi[self.rank]["src_global_x_overlap"].stop,
            self.map_mpi[self.rank]["src_global_y_overlap"].start,
            self.map_mpi[self.rank]["src_global_y_overlap"].stop)

        if self.horizontal_resampling:
            data = resample_2d_to_grid(self.read_axis_x(type="source_global", with_overlap=True),
                                       self.read_axis_y(type="source_global", with_overlap=True),
                                       self.read_axis_x(type="target", with_overlap=True),
                                       self.read_axis_y(type="target", with_overlap=True),
                                       data,
                                       Coverage.HORIZONTAL_INTERPOLATION_METHOD)

        return data[self.map_mpi[self.rank]["dst_local_y"], self.map_mpi[self.rank]["dst_local_x"]]

    def read_variable_sea_surface_wave_peak_period_at_time(self,t):

        index_t = self.find_time_index(t);

        data = self.reader.read_variable_sea_surface_wave_peak_period_at_time(
            self.map_mpi[self.rank]["src_global_t"].start + index_t,
            self.map_mpi[self.rank]["src_global_x_overlap"].start,
            self.map_mpi[self.rank]["src_global_x_overlap"].stop,
            self.map_mpi[self.rank]["src_global_y_overlap"].start,
            self.map_mpi[self.rank]["src_global_y_overlap"].stop)

        if self.horizontal_resampling:
            data = resample_2d_to_grid(self.read_axis_x(type="source_global", with_overlap=True),
                                       self.read_axis_y(type="source_global", with_overlap=True),
                                       self.read_axis_x(type="target", with_overlap=True),
                                       self.read_axis_y(type="target", with_overlap=True),
                                       data,
                                       Coverage.HORIZONTAL_INTERPOLATION_METHOD)

        return data[self.map_mpi[self.rank]["dst_local_y"], self.map_mpi[self.rank]["dst_local_x"]]

    def read_variable_sea_surface_wave_from_direction_at_time(self,t):

        index_t = self.find_time_index(t);

        data = self.reader.read_variable_sea_surface_wave_from_direction_at_time(
            self.map_mpi[self.rank]["src_global_t"].start + index_t,
            self.map_mpi[self.rank]["src_global_x_overlap"].start,
            self.map_mpi[self.rank]["src_global_x_overlap"].stop,
            self.map_mpi[self.rank]["src_global_y_overlap"].start,
            self.map_mpi[self.rank]["src_global_y_overlap"].stop)

        if self.horizontal_resampling:
            data = resample_2d_to_grid(self.read_axis_x(type="source_global", with_overlap=True),
                                       self.read_axis_y(type="source_global", with_overlap=True),
                                       self.read_axis_x(type="target", with_overlap=True),
                                       self.read_axis_y(type="target", with_overlap=True),
                                       data,
                                       Coverage.HORIZONTAL_INTERPOLATION_METHOD)

        return data[self.map_mpi[self.rank]["dst_local_y"], self.map_mpi[self.rank]["dst_local_x"]]

    def read_variable_sea_surface_wave_to_direction_at_time(self,t):

        index_t = self.find_time_index(t);

        data = self.reader.read_variable_sea_surface_wave_to_direction_at_time(
            self.map_mpi[self.rank]["src_global_t"].start + index_t,
            self.map_mpi[self.rank]["src_global_x_overlap"].start,
            self.map_mpi[self.rank]["src_global_x_overlap"].stop,
            self.map_mpi[self.rank]["src_global_y_overlap"].start,
            self.map_mpi[self.rank]["src_global_y_overlap"].stop)

        if self.horizontal_resampling:
            data = resample_2d_to_grid(self.read_axis_x(type="source_global", with_overlap=True),
                                       self.read_axis_y(type="source_global", with_overlap=True),
                                       self.read_axis_x(type="target", with_overlap=True),
                                       self.read_axis_y(type="target", with_overlap=True),
                                       data,
                                       Coverage.HORIZONTAL_INTERPOLATION_METHOD)

        return data[self.map_mpi[self.rank]["dst_local_y"], self.map_mpi[self.rank]["dst_local_x"]]

    def read_variable_sea_surface_wave_stokes_drift_velocity_at_time(self, t):
        """Retourne la dérive de Stokes en surface à la date souhaitée sur toute la couverture horizontale.
    @type t: datetime ou l'index
    @param t: date souhaitée
    @return: un tableau en deux dimensions [y,x]."""

        index_t = self.find_time_index(t);

        data = self.reader.read_variable_sea_surface_wave_stokes_drift_velocity_at_time(
            self.map_mpi[self.rank]["src_global_t"].start + index_t,
            self.map_mpi[self.rank]["src_global_x_overlap"].start,
            self.map_mpi[self.rank]["src_global_x_overlap"].stop,
            self.map_mpi[self.rank]["src_global_y_overlap"].start,
            self.map_mpi[self.rank]["src_global_y_overlap"].stop)

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

    def read_variable_radiation_pressure_bernouilli_head_at_time(self, t):
        """Retourne la pression J due aux vagues à la date souhaitée sur toute la couverture horizontale.
    @type t: datetime ou l'index
    @param t: date souhaitée
    @return: un tableau en deux dimensions [y,x]."""

        index_t = self.find_time_index(t);

        data = self.reader.read_variable_radiation_pressure_bernouilli_head_at_time(
            self.map_mpi[self.rank]["src_global_t"].start + index_t,
            self.map_mpi[self.rank]["src_global_x_overlap"].start,
            self.map_mpi[self.rank]["src_global_x_overlap"].stop,
            self.map_mpi[self.rank]["src_global_y_overlap"].start,
            self.map_mpi[self.rank]["src_global_y_overlap"].stop)

        if self.horizontal_resampling:
            data = resample_2d_to_grid(self.read_axis_x(type="source_global", with_overlap=True),
                                       self.read_axis_y(type="source_global", with_overlap=True),
                                       self.read_axis_x(type="target", with_overlap=True),
                                       self.read_axis_y(type="target", with_overlap=True),
                                       data,
                                       Coverage.HORIZONTAL_INTERPOLATION_METHOD)

        return data[self.map_mpi[self.rank]["dst_local_y"], self.map_mpi[self.rank]["dst_local_x"]]

    def read_variable_sea_surface_wave_energy_flux_to_ocean_at_time(self, t):
        """Retourne la waves_to_ocean_energy_flux à la date souhaitée sur toute la couverture horizontale.
    @type t: datetime ou l'index
    @param t: date souhaitée
    @return: un tableau en deux dimensions [y,x]."""

        index_t = self.find_time_index(t);

        data = self.reader.read_variable_sea_surface_wave_energy_flux_to_ocean_at_time(
            self.map_mpi[self.rank]["src_global_t"].start + index_t,
            self.map_mpi[self.rank]["src_global_x_overlap"].start,
            self.map_mpi[self.rank]["src_global_x_overlap"].stop,
            self.map_mpi[self.rank]["src_global_y_overlap"].start,
            self.map_mpi[self.rank]["src_global_y_overlap"].stop)

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

    def read_variable_sea_surface_wave_energy_dissipation_at_ground_level_at_time(self, t):
        """Retourne la l'énergie des vagues dissipée par le fond à la date souhaitée sur toute la couverture horizontale.
    @type t: datetime ou l'index
    @param t: date souhaitée
    @return: un tableau en deux dimensions [y,x]."""

        index_t = self.find_time_index(t);

        data = self.reader.read_variable_sea_surface_wave_energy_dissipation_at_ground_level_at_time(
            self.map_mpi[self.rank]["src_global_t"].start + index_t,
            self.map_mpi[self.rank]["src_global_x_overlap"].start,
            self.map_mpi[self.rank]["src_global_x_overlap"].stop,
            self.map_mpi[self.rank]["src_global_y_overlap"].start,
            self.map_mpi[self.rank]["src_global_y_overlap"].stop)

        if self.horizontal_resampling:
            data = resample_2d_to_grid(self.read_axis_x(type="source_global", with_overlap=True),
                                       self.read_axis_y(type="source_global", with_overlap=True),
                                       self.read_axis_x(type="target", with_overlap=True),
                                       self.read_axis_y(type="target", with_overlap=True),
                                       data,
                                       Coverage.HORIZONTAL_INTERPOLATION_METHOD)

        return data[self.map_mpi[self.rank]["dst_local_y"], self.map_mpi[self.rank]["dst_local_x"]]

    #################
    # WAVES
    # Momentum flux
    #################
    def read_variable_atmosphere_momentum_flux_to_waves_at_time(self,t):
        """Retourne la composante u du tau atmosphere->vagues à la date souhaitée sur toute la couverture horizontale.
    @type t: datetime ou l'index
    @param t: date souhaitée
    @return: un tableau en deux dimensions [y,x]."""

        index_t = self.find_time_index(t);

        data = self.reader.read_variable_atmosphere_momentum_flux_to_waves_at_time(
            self.map_mpi[self.rank]["src_global_t"].start + index_t,
            self.map_mpi[self.rank]["src_global_x_overlap"].start,
            self.map_mpi[self.rank]["src_global_x_overlap"].stop,
            self.map_mpi[self.rank]["src_global_y_overlap"].start,
            self.map_mpi[self.rank]["src_global_y_overlap"].stop)

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

    def read_variable_waves_momentum_flux_to_ocean_at_time(self,t):
        """Retourne la composante u du tau vagues->ocean à la date souhaitée sur toute la couverture horizontale.
    @type t: datetime ou l'index
    @param t: date souhaitée
    @return: un tableau en deux dimensions [y,x]."""

        index_t = self.find_time_index(t);

        data = self.reader.read_variable_waves_momentum_flux_to_ocean_at_time(
            self.map_mpi[self.rank]["src_global_t"].start + index_t,
            self.map_mpi[self.rank]["src_global_x_overlap"].start,
            self.map_mpi[self.rank]["src_global_x_overlap"].stop,
            self.map_mpi[self.rank]["src_global_y_overlap"].start,
            self.map_mpi[self.rank]["src_global_y_overlap"].stop)

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

    #################
    # METEO
    # 2D
    #################
    def read_variable_rainfall_amount_at_time(self, t):
        """Retourne les composantes u,v de rain à la date souhaitée
    @type t: datetime ou l'index
    @param t: date souhaitée
    @return: un tableau en deux dimensions [u_comp,v_comp] contenant chacun deux dimensions [y,x]."""

        index_t = self.find_time_index(t);

        data = self.reader.read_variable_rainfall_amount_at_time(
            self.map_mpi[self.rank]["src_global_t"].start + index_t,
            self.map_mpi[self.rank]["src_global_x_overlap"].start,
            self.map_mpi[self.rank]["src_global_x_overlap"].stop,
            self.map_mpi[self.rank]["src_global_y_overlap"].start,
            self.map_mpi[self.rank]["src_global_y_overlap"].stop)

        if self.horizontal_resampling:
            data = resample_2d_to_grid(self.read_axis_x(type="source_global", with_overlap=True),
                                       self.read_axis_y(type="source_global", with_overlap=True),
                                       self.read_axis_x(type="target", with_overlap=True),
                                       self.read_axis_y(type="target", with_overlap=True),
                                       data,
                                       Coverage.HORIZONTAL_INTERPOLATION_METHOD)

        return data[self.map_mpi[self.rank]["dst_local_y"], self.map_mpi[self.rank]["dst_local_x"]]

    #################
    # METEO
    # Sea surface
    #################

    def read_variable_surface_air_pressure_at_time(self,t):
        """Retourne la pression à la surface à la date souhaitée sur toute la couverture horizontale.
    @type t: datetime ou l'index
    @param t: date souhaitée
    @return: un tableau en deux dimensions [y,x]."""

        index_t = self.find_time_index(t);

        data = self.reader.read_variable_surface_air_pressure_at_time(
            self.map_mpi[self.rank]["src_global_t"].start + index_t,
            self.map_mpi[self.rank]["src_global_x_overlap"].start,
            self.map_mpi[self.rank]["src_global_x_overlap"].stop,
            self.map_mpi[self.rank]["src_global_y_overlap"].start,
            self.map_mpi[self.rank]["src_global_y_overlap"].stop)

        if self.horizontal_resampling:
            data = resample_2d_to_grid(self.read_axis_x(type="source_global", with_overlap=True),
                                       self.read_axis_y(type="source_global", with_overlap=True),
                                       self.read_axis_x(type="target", with_overlap=True),
                                       self.read_axis_y(type="target", with_overlap=True),
                                       data,
                                       Coverage.HORIZONTAL_INTERPOLATION_METHOD)

        return data[self.map_mpi[self.rank]["dst_local_y"], self.map_mpi[self.rank]["dst_local_x"]]

    def read_variable_sea_surface_air_pressure_at_time(self, t):
        """Retourne la pression à la surface à la date souhaitée sur toute la couverture horizontale.
    @type t: datetime ou l'index
    @param t: date souhaitée
    @return: un tableau en deux dimensions [y,x]."""

        index_t = self.find_time_index(t);

        data = self.reader.read_variable_sea_surface_air_pressure_at_time(
            self.map_mpi[self.rank]["src_global_t"].start + index_t,
            self.map_mpi[self.rank]["src_global_x_overlap"].start,
            self.map_mpi[self.rank]["src_global_x_overlap"].stop,
            self.map_mpi[self.rank]["src_global_y_overlap"].start,
            self.map_mpi[self.rank]["src_global_y_overlap"].stop)

        if self.horizontal_resampling:
            data = resample_2d_to_grid(self.read_axis_x(type="source_global", with_overlap=True),
                                       self.read_axis_y(type="source_global", with_overlap=True),
                                       self.read_axis_x(type="target", with_overlap=True),
                                       self.read_axis_y(type="target", with_overlap=True),
                                       data,
                                       Coverage.HORIZONTAL_INTERPOLATION_METHOD)

        return data[self.map_mpi[self.rank]["dst_local_y"], self.map_mpi[self.rank]["dst_local_x"]]

    def read_variable_wind_stress_at_time(self, t):
        """Retourne les composantes u,v de la contrainte de vent à la date souhaitée
    @type t: datetime ou l'index
    @param t: date souhaitée
    @return: un tableau en deux dimensions [u_comp,v_comp] contenant chacun deux dimensions [y,x]."""

        index_t = self.find_time_index(t);

        data = self.reader.read_variable_wind_stress_at_time(
            self.map_mpi[self.rank]["src_global_t"].start + index_t,
            self.map_mpi[self.rank]["src_global_x_overlap"].start,
            self.map_mpi[self.rank]["src_global_x_overlap"].stop,
            self.map_mpi[self.rank]["src_global_y_overlap"].start,
            self.map_mpi[self.rank]["src_global_y_overlap"].stop)

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

    def read_variable_surface_downward_sensible_heat_flux_at_time(self, t):
        """Retourne les composantes u,v de surface sensible heat flux à la date souhaitée
    @type t: datetime ou l'index
    @param t: date souhaitée
    @return: un tableau en deux dimensions [u_comp,v_comp] contenant chacun deux dimensions [y,x]."""

        index_t = self.find_time_index(t);

        data = self.rreader.read_variable_surface_downward_sensible_heat_flux_at_time(
            self.map_mpi[self.rank]["src_global_t"].start + index_t,
            self.map_mpi[self.rank]["src_global_x_overlap"].start,
            self.map_mpi[self.rank]["src_global_x_overlap"].stop,
            self.map_mpi[self.rank]["src_global_y_overlap"].start,
            self.map_mpi[self.rank]["src_global_y_overlap"].stop)

        if self.horizontal_resampling:
            data = resample_2d_to_grid(self.read_axis_x(type="source_global", with_overlap=True),
                                       self.read_axis_y(type="source_global", with_overlap=True),
                                       self.read_axis_x(type="target", with_overlap=True),
                                       self.read_axis_y(type="target", with_overlap=True),
                                       data,
                                       Coverage.HORIZONTAL_INTERPOLATION_METHOD)

        return data[self.map_mpi[self.rank]["dst_local_y"], self.map_mpi[self.rank]["dst_local_x"]]

    def read_variable_surface_downward_latent_heat_flux_at_time(self, t):
        """Retourne les composantes u,v de surface latente heat flux à la date souhaitée
    @type t: datetime ou l'index
    @param t: date souhaitée
    @return: un tableau en deux dimensions [u_comp,v_comp] contenant chacun deux dimensions [y,x]."""

        index_t = self.find_time_index(t);

        data = self.reader.read_variable_surface_downward_latent_heat_flux_at_time(
            self.map_mpi[self.rank]["src_global_t"].start + index_t,
            self.map_mpi[self.rank]["src_global_x_overlap"].start,
            self.map_mpi[self.rank]["src_global_x_overlap"].stop,
            self.map_mpi[self.rank]["src_global_y_overlap"].start,
            self.map_mpi[self.rank]["src_global_y_overlap"].stop)

        if self.horizontal_resampling:
            data = resample_2d_to_grid(self.read_axis_x(type="source_global", with_overlap=True),
                                       self.read_axis_y(type="source_global", with_overlap=True),
                                       self.read_axis_x(type="target", with_overlap=True),
                                       self.read_axis_y(type="target", with_overlap=True),
                                       data,
                                       Coverage.HORIZONTAL_INTERPOLATION_METHOD)

        return data[self.map_mpi[self.rank]["dst_local_y"], self.map_mpi[self.rank]["dst_local_x"]]

    def read_variable_surface_air_temperature_at_time(self, t):
        """Retourne les composantes u,v de surface air temperature à la date souhaitée
    @type t: datetime ou l'index
    @param t: date souhaitée
    @return: un tableau en deux dimensions [u_comp,v_comp] contenant chacun deux dimensions [y,x]."""

        index_t = self.find_time_index(t);

        data = self.reader.read_variable_surface_air_temperature_at_time(
            self.map_mpi[self.rank]["src_global_t"].start + index_t,
            self.map_mpi[self.rank]["src_global_x_overlap"].start,
            self.map_mpi[self.rank]["src_global_x_overlap"].stop,
            self.map_mpi[self.rank]["src_global_y_overlap"].start,
            self.map_mpi[self.rank]["src_global_y_overlap"].stop)

        if self.horizontal_resampling:
            data = resample_2d_to_grid(self.read_axis_x(type="source_global", with_overlap=True),
                                       self.read_axis_y(type="source_global", with_overlap=True),
                                       self.read_axis_x(type="target", with_overlap=True),
                                       self.read_axis_y(type="target", with_overlap=True),
                                       data,
                                       Coverage.HORIZONTAL_INTERPOLATION_METHOD)

        return data[self.map_mpi[self.rank]["dst_local_y"], self.map_mpi[self.rank]["dst_local_x"]]

    def read_variable_dew_point_temperature_at_time(self, t):
        """Retourne les composantes u,v de dewpoint temperature à la date souhaitée
    @type t: datetime ou l'index
    @param t: date souhaitée
    @return: un tableau en deux dimensions [u_comp,v_comp] contenant chacun deux dimensions [y,x]."""

        index_t = self.find_time_index(t);

        data = self.reader.read_variable_dew_point_temperature_at_time(
            self.map_mpi[self.rank]["src_global_t"].start + index_t,
            self.map_mpi[self.rank]["src_global_x_overlap"].start,
            self.map_mpi[self.rank]["src_global_x_overlap"].stop,
            self.map_mpi[self.rank]["src_global_y_overlap"].start,
            self.map_mpi[self.rank]["src_global_y_overlap"].stop)

        if self.horizontal_resampling:
            data = resample_2d_to_grid(self.read_axis_x(type="source_global", with_overlap=True),
                                       self.read_axis_y(type="source_global", with_overlap=True),
                                       self.read_axis_x(type="target", with_overlap=True),
                                       self.read_axis_y(type="target", with_overlap=True),
                                       data,
                                       Coverage.HORIZONTAL_INTERPOLATION_METHOD)

        return data[self.map_mpi[self.rank]["dst_local_y"], self.map_mpi[self.rank]["dst_local_x"]]

    def read_variable_surface_downwards_solar_radiation_at_time(self, t):
        """Retourne les composantes u,v de surface solar radiation downwards à la date souhaitée
    @type t: datetime ou l'index
    @param t: date souhaitée
    @return: un tableau en deux dimensions [u_comp,v_comp] contenant chacun deux dimensions [y,x]."""

        index_t = self.find_time_index(t);

        data = self.reader.read_variable_surface_downwards_solar_radiation_at_time(
            self.map_mpi[self.rank]["src_global_t"].start + index_t,
            self.map_mpi[self.rank]["src_global_x_overlap"].start,
            self.map_mpi[self.rank]["src_global_x_overlap"].stop,
            self.map_mpi[self.rank]["src_global_y_overlap"].start,
            self.map_mpi[self.rank]["src_global_y_overlap"].stop)

        if self.horizontal_resampling:
            data = resample_2d_to_grid(self.read_axis_x(type="source_global", with_overlap=True),
                                       self.read_axis_y(type="source_global", with_overlap=True),
                                       self.read_axis_x(type="target", with_overlap=True),
                                       self.read_axis_y(type="target", with_overlap=True),
                                       data,
                                       Coverage.HORIZONTAL_INTERPOLATION_METHOD)

        return data[self.map_mpi[self.rank]["dst_local_y"], self.map_mpi[self.rank]["dst_local_x"]]

    def read_variable_surface_downwards_thermal_radiation_at_time(self, t):
        """Retourne les composantes u,v de surface thermal radiation downwards à la date souhaitée
    @type t: datetime ou l'index
    @param t: date souhaitée
    @return: un tableau en deux dimensions [u_comp,v_comp] contenant chacun deux dimensions [y,x]."""

        index_t = self.find_time_index(t);

        data = self.reader.read_variable_surface_downwards_thermal_radiation_at_time(
            self.map_mpi[self.rank]["src_global_t"].start + index_t,
            self.map_mpi[self.rank]["src_global_x_overlap"].start,
            self.map_mpi[self.rank]["src_global_x_overlap"].stop,
            self.map_mpi[self.rank]["src_global_y_overlap"].start,
            self.map_mpi[self.rank]["src_global_y_overlap"].stop)

        if self.horizontal_resampling:
            data = resample_2d_to_grid(self.read_axis_x(type="source_global", with_overlap=True),
                                       self.read_axis_y(type="source_global", with_overlap=True),
                                       self.read_axis_x(type="target", with_overlap=True),
                                       self.read_axis_y(type="target", with_overlap=True),
                                       data,
                                       Coverage.HORIZONTAL_INTERPOLATION_METHOD)

        return data[self.map_mpi[self.rank]["dst_local_y"], self.map_mpi[self.rank]["dst_local_x"]]

    def read_variable_surface_solar_radiation_at_time(self, t):
        """Retourne les composantes u,v de surface solar radiation à la date souhaitée
    @type t: datetime ou l'index
    @param t: date souhaitée
    @return: un tableau en deux dimensions [u_comp,v_comp] contenant chacun deux dimensions [y,x]."""

        index_t = self.find_time_index(t);

        data = self.reader.read_variable_surface_solar_radiation_at_time(
            self.map_mpi[self.rank]["src_global_t"].start + index_t,
            self.map_mpi[self.rank]["src_global_x_overlap"].start,
            self.map_mpi[self.rank]["src_global_x_overlap"].stop,
            self.map_mpi[self.rank]["src_global_y_overlap"].start,
            self.map_mpi[self.rank]["src_global_y_overlap"].stop)

        if self.horizontal_resampling:
            data = resample_2d_to_grid(self.read_axis_x(type="source_global", with_overlap=True),
                                       self.read_axis_y(type="source_global", with_overlap=True),
                                       self.read_axis_x(type="target", with_overlap=True),
                                       self.read_axis_y(type="target", with_overlap=True),
                                       data,
                                       Coverage.HORIZONTAL_INTERPOLATION_METHOD)

        return data[self.map_mpi[self.rank]["dst_local_y"], self.map_mpi[self.rank]["dst_local_x"]]

    def read_variable_surface_thermal_radiation_at_time(self, t):
        """Retourne les composantes u,v de surface thermal radiation à la date souhaitée
    @type t: datetime ou l'index
    @param t: date souhaitée
    @return: un tableau en deux dimensions [u_comp,v_comp] contenant chacun deux dimensions [y,x]."""

        index_t = self.find_time_index(t);

        data = self.reader.read_variable_surface_thermal_radiation_at_time(
            self.map_mpi[self.rank]["src_global_t"].start + index_t,
            self.map_mpi[self.rank]["src_global_x_overlap"].start,
            self.map_mpi[self.rank]["src_global_x_overlap"].stop,
            self.map_mpi[self.rank]["src_global_y_overlap"].start,
            self.map_mpi[self.rank]["src_global_y_overlap"].stop)

        if self.horizontal_resampling:
            data = resample_2d_to_grid(self.read_axis_x(type="source_global", with_overlap=True),
                                       self.read_axis_y(type="source_global", with_overlap=True),
                                       self.read_axis_x(type="target", with_overlap=True),
                                       self.read_axis_y(type="target", with_overlap=True),
                                       data,
                                       Coverage.HORIZONTAL_INTERPOLATION_METHOD)

        return data[self.map_mpi[self.rank]["dst_local_y"], self.map_mpi[self.rank]["dst_local_x"]]

    #################
    # METEO
    # At 10 m
    #################
    def read_variable_wind_10m_at_time(self,t):
        """Retourne les composantes u,v du vent à la date souhaitée
    @type t: datetime ou l'index
    @param t: date souhaitée
    @return: un tableau en deux dimensions [u_comp,v_comp] contenant chacun deux dimensions [y,x]."""

        index_t = self.find_time_index(t);

        data = self.reader.read_variable_wind_10m_at_time(
            self.map_mpi[self.rank]["src_global_t"].start + index_t,
            self.map_mpi[self.rank]["src_global_x_overlap"].start,
            self.map_mpi[self.rank]["src_global_x_overlap"].stop,
            self.map_mpi[self.rank]["src_global_y_overlap"].start,
            self.map_mpi[self.rank]["src_global_y_overlap"].stop)

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

    def read_variable_wind_speed_10m_at_time(self, date):
        comp = self.read_variable_wind_10m_at_time(date)
        result = np.zeros([self.get_y_size(),self.get_x_size()])
        result[:] = np.nan
        for x in range(0,self.get_x_size()):
            for y in range(0,self.get_y_size()):
                result[y,x] = math.sqrt(comp[0][y,x] ** 2 + comp[1][y,x] ** 2)

        return result

    def read_variable_wind_from_direction_10m_at_time(self, date):
        comp = self.read_variable_wind_10m_at_time(date)
        result = np.zeros([self.get_y_size(), self.get_x_size()])
        result[:] = np.nan
        for x in range(0, self.get_x_size()):
            for y in range(0, self.get_y_size()):
                result[y, x] = 270. - (180.0 / math.pi) * (math.atan2(comp[0][y,x], comp[1][y,x])) + 180.0 % 360.0

        return result

    def read_variable_wind_to_direction_10m_at_time(self, date):
        comp = self.read_variable_wind_10m_at_time(date)
        result = np.zeros([self.get_y_size(), self.get_x_size()])
        result[:] = np.nan
        for x in range(0, self.get_x_size()):
            for y in range(0, self.get_y_size()):
                result[y, x] = 270. - (180.0 / math.pi) * (math.atan2(comp[0][y,x], comp[1][y,x])) % 360.0

        return result









