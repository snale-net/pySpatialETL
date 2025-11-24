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
import math
import os
from concurrent.futures import ProcessPoolExecutor
from datetime import datetime
from datetime import timedelta

import cftime
import numpy as np
import pandas
from array_split import shape_split

from spatialetl.coverage.coverage import Coverage
from spatialetl.exception.not_found_in_rank_error import NotFoundInRankError
from spatialetl.operator.interpolator.interpolator_core import resample_2d_to_grid, resample_faster_2d_to_grid
from spatialetl.utils.logger import logging


class TimeCoverage(Coverage):
    """
    La classe TimeCoverage est une extension de la classe Coverage.
    Elle rajoute une dimension temporelle à la couverture horizontale classique.
    """

    TIME_DATUM = datetime(1970, 1, 1)
    TIME_DELTA = timedelta(minutes=15)
    TIME_OVERLAPING_SIZE = 0

    def __init__(self, reader, bbox=None, resolution_x=None, resolution_y=None, start_time=None, end_time=None,
                 freq=None,nb_thread:int=os.cpu_count()-1):

        Coverage.__init__(self, reader, bbox=bbox, resolution_x=resolution_x, resolution_y=resolution_y, nb_thread=nb_thread);

        self.source_global_t_size = self.reader.get_t_size()
        self.source_global_axis_t = self.reader.read_axis_t(0, self.source_global_t_size, 0);

        self.temporal_resampling = False
        tmin = 0
        tmax = self.source_global_t_size
        zero_delta = timedelta(minutes=00)

        if start_time is not None:

            if type(start_time) == datetime:
                time = start_time
            elif type(start_time) == str:
                try:
                    time = datetime.strptime(start_time, '%Y-%m-%d %H:%M:%S')
                except ValueError as ex:
                    raise ValueError("start_time is not well formated " + str(ex))
            else:
                raise ValueError("start_time have to be string or datetime. Found " + str(type(start_time)))

            nearest_t_index = (np.abs(np.asarray(self.source_global_axis_t) - time)).argmin()

            if time - datetime.strptime(str(self.source_global_axis_t[nearest_t_index]),
                                        '%Y-%m-%d %H:%M:%S') == zero_delta or abs(
                    time - datetime.strptime(str(self.source_global_axis_t[nearest_t_index]),
                                             '%Y-%m-%d %H:%M:%S')) < TimeCoverage.TIME_DELTA:
                tmin = nearest_t_index
            else:
                raise ValueError(str(time) + " not found. Maybe the TimeCoverage.TIME_DELTA (" + str(
                    TimeCoverage.TIME_DELTA) + ") is too small or the date is out the range.")

        if end_time is not None:

            if type(end_time) == datetime:
                time = end_time
            elif type(end_time) == str:
                try:
                    time = datetime.strptime(end_time, '%Y-%m-%d %H:%M:%S')
                except ValueError as ex:
                    raise ValueError("end_time is not well formated " + str(ex))
            else:
                raise ValueError("end_time have to be string or datetime. Found " + str(type(end_time)))

            nearest_t_index = (np.abs(np.asarray(self.source_global_axis_t) - time)).argmin()

            if time - datetime.strptime(str(self.source_global_axis_t[nearest_t_index]),
                                        '%Y-%m-%d %H:%M:%S') == zero_delta or abs(
                    time - datetime.strptime(str(self.source_global_axis_t[nearest_t_index]),
                                             '%Y-%m-%d %H:%M:%S')) < TimeCoverage.TIME_DELTA:
                tmax = nearest_t_index + 1
            else:
                raise ValueError(str(time) + " not found. Maybe the TimeCoverage.TIME_DELTA (" + str(
                    TimeCoverage.TIME_DELTA) + ") is too small or the date is out the range.")

        if freq is not None:
            self.target_global_axis_t = pandas.date_range(start=datetime.utcfromtimestamp(
                self.read_axis_t(type="source_global", with_overlap=False, timestamp=1)[tmin]),
                                                          end=datetime.utcfromtimestamp(
                                                              self.read_axis_t(type="source_global", with_overlap=False,
                                                                               timestamp=1)[tmax - 1]),
                                                          freq=freq).to_pydatetime();
            self.target_global_t_size = np.shape(self.target_global_axis_t)[0]
        else:
            self.target_global_axis_t = self.source_global_axis_t[tmin:tmax]
            self.target_global_t_size = tmax - tmin

        self.__init_parallel_map()

        if type(self) == TimeCoverage and self.horizontal_resampling and self.rank == 0:
            logging.info(
                '[horizontal_interpolation] Source grid size : (' + str(self.source_global_x_size) + ", " + str(
                    self.source_global_y_size) + ")")
            logging.info(
                '[horizontal_interpolation] Target grid size : (' + str(self.target_global_x_size) + ", " + str(
                    self.target_global_y_size) + ")")

            Coverage.compute_weight(self)

        if self.rank == 0 and self.comm:
            logging.debug("MPI map:")

        if self.comm:
            logging.debug(f"{"-" * 10} Proc n° {self.rank} {"-" * 10}")
        else:
            logging.debug("Multithreads map:")
            logging.debug(f"{"-" * 10} Target grid {"-" * 10}")
            for key in ['dst_global_t','dst_global_x', 'dst_global_y', 'dst_local_t_size', 'dst_local_x_size', 'dst_local_y_size']:
                logging.debug(f"{key} = {self.parallel_map[self.rank][key]}")

        for key in self.parallel_map[self.rank]:
            if self.comm:
                logging.debug(f"    {key} = {self.parallel_map[self.rank][key]}")

        for thread in range(len(self.parallel_map[self.rank]["threads_map"])):
            logging.debug(f"   {"-" * 10} Thread n° {thread} {"-" * 10}")
            for thread_key in self.parallel_map[self.rank]["threads_map"][thread]:
                logging.debug(f"    {thread_key} = {self.parallel_map[self.rank]["threads_map"][thread][thread_key]}")

        if self.rank == 0:
            logging.debug("-" * 20)

    def __init_parallel_map(self):
        """
        Create the MPI map for parallel processing.
        The MPI map is a dictionary that contains the mapping of the source and destination grids for each MPI rank.

        Examples
        --------
        >>> coverage.__init_parellel_map()
        """

        # Découpage des axes
        # if self.horizontal_resampling:
        #     # Découpage sur le temps uniquemenent
        #     target_slices = shape_split(target_sample, self.size, axis=[0, 1, 1])
        #     # Si on feet pas le nombre de proc
        #     if len( target_slices.flatten()) != self.size:
        #         target_slices = shape_split(target_sample, self.size, axis=[0, 0, 0])
        # else:

        target_mpi_sample = (self.target_global_t_size, self.target_global_y_size, self.target_global_x_size)

        # Split the axes with the MPI size
        target_mpi_slices = shape_split(target_mpi_sample, self.size, axis=[0, 0, 0])

        mpi_slice_index = 0
        for slyce in target_mpi_slices.flatten():
            mpi_slice = tuple(slyce)
            self.parallel_map[mpi_slice_index] = self.__compute_mpi_slice(mpi_slice)

            # Split the axes with the number of threads
            target_threads_sample = (self.parallel_map[mpi_slice_index]["dst_local_t_size"],self.parallel_map[mpi_slice_index]["dst_local_y_size"],
                                     self.parallel_map[mpi_slice_index]["dst_local_x_size"])
            target_threads_slices = shape_split(target_threads_sample, self.threads_number, axis=[0, 0, 0])

            # If we can't divide the grid dimensions with the number of threads,
            # we set the threads number with the number of slices
            self.threads_number = len(target_threads_slices.flatten())

            self.parallel_map[mpi_slice_index]['threads_map'] = np.empty([self.threads_number], dtype=object)

            slice_thread_index = 0
            for thread_slyce in target_threads_slices.flatten():
                thread_slice = tuple(thread_slyce)
                self.parallel_map[mpi_slice_index]["threads_map"][slice_thread_index] = self.__compute_thread_slice(
                    thread_slice)
                slice_thread_index = slice_thread_index + 1

            mpi_slice_index = mpi_slice_index + 1

        self.__update_mpi_map()
        self.__update_thread_map()

    def __compute_mpi_slice(self,slice):

        map = Coverage.compute_mpi_slice(self,slice[1:])

        # Grille source
        map["dst_global_t"] = slice[0]

        map["dst_local_t_size"] = map["dst_global_t"].stop - map["dst_global_t"].start

        dst_global_t_min_overlap = max(0, map["dst_global_t"].start - TimeCoverage.TIME_OVERLAPING_SIZE)
        dst_global_t_max_overlap = min(self.target_global_t_size,
                                       map["dst_global_t"].stop + TimeCoverage.TIME_OVERLAPING_SIZE)
        map["dst_global_t_overlap"] = np.s_[dst_global_t_min_overlap:dst_global_t_max_overlap]

        map["dst_global_t_size_overlap"] = map["dst_global_t_overlap"].stop - map["dst_global_t_overlap"].start

        dst_t_min = TimeCoverage.TIME_OVERLAPING_SIZE
        dst_t_max = map["dst_global_t_size_overlap"] - TimeCoverage.TIME_OVERLAPING_SIZE

        if map["dst_global_t"].start == 0:
            dst_t_min = 0

        if map["dst_global_t"].stop == self.target_global_t_size:
            dst_t_max = map["dst_global_t_size_overlap"]

        map["dst_local_t"] = np.s_[dst_t_min:dst_t_max]

        # Source grille
        map["src_global_t"] = map["dst_global_t"]

        map["src_global_t_overlap"] = map["dst_global_t_overlap"]

        map["src_local_t"] = map["dst_local_t"]

        map["src_local_t_size"] = map["dst_local_t_size"]

        map["src_local_t_size_overlap"] = map["dst_global_t_size_overlap"]

        return map

    def __compute_thread_slice(self, slice):
        """
        Create the multi threading map for parallel processing.

        Examples
        --------
        >>> coverage.__create_threading_map()
        """
        map = Coverage.compute_thread_slice(self,slice[1:])

        #### Destination grid ###
        map["dst_global_t"] = slice[0]

        map["dst_local_t_size"] = map["dst_global_t"].stop - map["dst_global_t"].start

        # Compute overlap
        dst_global_t_min_overlap = max(0, map["dst_global_t"].start - TimeCoverage.TIME_OVERLAPING_SIZE)
        dst_global_t_max_overlap = min(self.parallel_map[self.rank]["dst_local_t_size"],
                                       map["dst_global_t"].stop + TimeCoverage.TIME_OVERLAPING_SIZE)
        map["dst_global_t_overlap"] = np.s_[dst_global_t_min_overlap:dst_global_t_max_overlap]


        map["dst_global_t_size_overlap"] = map["dst_global_t_overlap"].stop - map["dst_global_t_overlap"].start

        # Compute dest local grid
        if map["dst_global_t"].start == 0:
            dst_local_t_min = 0
        elif map["dst_global_t"].start == int(TimeCoverage.TIME_OVERLAPING_SIZE / 2):
            dst_local_t_min = int(TimeCoverage.TIME_OVERLAPING_SIZE / 2)
        elif map["dst_global_t"].start == TimeCoverage.TIME_OVERLAPING_SIZE:
            dst_local_t_min = TimeCoverage.TIME_OVERLAPING_SIZE
        else:
            dst_local_t_min = map["dst_global_t"].start - map["dst_global_t_overlap"].start

        dst_local_t_max = dst_local_t_min + map["dst_local_t_size"]

        map["dst_local_t"] = np.s_[dst_local_t_min:dst_local_t_max]

        ### Source grille ###
        map["src_global_t"] = map["dst_global_t"]
        map["src_global_t_overlap"] = map["dst_global_t_overlap"]
        map["src_local_t"] = map["dst_local_t"]
        map["src_local_t_size"] = map["dst_local_t_size"]
        map["src_local_t_size_overlap"] = map["dst_global_t_size_overlap"]

        return map

    def __update_mpi_map(self):

        Coverage.update_mpi_map(self)

        if self.get_t_size(type="target", with_overlap=False) == 1:
            tmin = (np.abs(np.asarray(self.read_axis_t(type="source_global", with_overlap=False, timestamp=1)) - np.min(
                self.read_axis_t(type="target", with_overlap=False, timestamp=1)))).argmin()
            tmax = tmin + 1
        else:
            idx = np.where((self.read_axis_t(type="source_global", with_overlap=False, timestamp=1) >= np.min(
                self.read_axis_t(type="target", with_overlap=False, timestamp=1))) &
                           (self.read_axis_t(type="source_global", with_overlap=False, timestamp=1) <= np.max(
                               self.read_axis_t(type="target", with_overlap=False, timestamp=1))))

            tmin = np.min(idx[0])
            tmax = np.max(idx[0]) + 1

        # SRC GLOBAL
        self.parallel_map[self.rank]["src_global_t"] = np.s_[tmin:tmax]
        self.parallel_map[self.rank]["src_global_t_size"] = tmax - tmin

        dst_global_t_min_overlap = max(0, self.parallel_map[self.rank][
            "src_global_t"].start - TimeCoverage.TIME_OVERLAPING_SIZE)
        dst_global_t_max_overlap = min(self.source_global_t_size,
                                       self.parallel_map[self.rank][
                                           "src_global_t"].stop + TimeCoverage.TIME_OVERLAPING_SIZE)
        self.parallel_map[self.rank]["src_global_t_overlap"] = np.s_[
                                                          dst_global_t_min_overlap:dst_global_t_max_overlap]

        self.parallel_map[self.rank]["src_global_t_size_overlap"] = self.parallel_map[self.rank][
                                                                   "src_global_t_overlap"].stop - \
                                                               self.parallel_map[self.rank][
                                                                   "src_global_t_overlap"].start

        self.parallel_map[self.rank]["src_local_t_size"] = tmax - tmin
        self.parallel_map[self.rank]["src_local_t"] = np.s_[0:self.parallel_map[self.rank]["src_local_t_size"]]

        # OVERLAP
        self.parallel_map[self.rank]["src_local_t_size_overlap"] = self.parallel_map[self.rank][
            "src_global_t_size_overlap"]

        self.parallel_map[self.rank]["src_local_t_overlap"] = np.s_[
                                                         0:self.parallel_map[self.rank]["src_local_t_size_overlap"]]

    def __update_thread_map(self):

        Coverage.update_thread_map(self)

        for thread in range(0, self.threads_number):

            if self.get_t_size(type="target", with_overlap=False) == 1:
                tmin = (np.abs(np.asarray(self.source_global_axis_t[self.parallel_map[self.rank]["src_global_t"]]) - np.min(
                    self.__read_thread_axis_t(type="target", with_overlap=False, timestamp=0, current_thread=thread)))).argmin()
                tmax = tmin + 1
            else:
                idx = np.where((np.asarray(self.source_global_axis_t[self.parallel_map[self.rank]["src_global_t"]]) >= np.min(
                    self.__read_thread_axis_t(type="target", with_overlap=False, timestamp=0, current_thread=thread))) &
                               (np.asarray(self.source_global_axis_t[self.parallel_map[self.rank]["src_global_t"]]) <= np.max(
                                   self.__read_thread_axis_t(type="target", with_overlap=False, timestamp=0, current_thread=thread))))

                tmin = np.min(idx[0])
                tmax = np.max(idx[0]) + 1

            # SRC GLOBAL
            self.parallel_map[self.rank]["threads_map"][thread]["src_global_t"] = np.s_[tmin:tmax]
            self.parallel_map[self.rank]["threads_map"][thread]["src_global_t_size"] = tmax - tmin

            dst_global_t_min_overlap = max(0, self.parallel_map[self.rank]["threads_map"][thread][
                "src_global_t"].start - TimeCoverage.TIME_OVERLAPING_SIZE)
            dst_global_t_max_overlap = min(self.source_global_t_size,
                                           self.parallel_map[self.rank]["threads_map"][thread][
                                               "src_global_t"].stop + TimeCoverage.TIME_OVERLAPING_SIZE)
            self.parallel_map[self.rank]["threads_map"][thread]["src_global_t_overlap"] = np.s_[
                dst_global_t_min_overlap:dst_global_t_max_overlap]

            self.parallel_map[self.rank]["threads_map"][thread]["src_global_t_size_overlap"] = self.parallel_map[self.rank]["threads_map"][thread][
                                                                            "src_global_t_overlap"].stop - \
                                                                        self.parallel_map[self.rank]["threads_map"][thread][
                                                                            "src_global_t_overlap"].start

            self.parallel_map[self.rank]["threads_map"][thread]["src_local_t_size"] = tmax - tmin
            self.parallel_map[self.rank]["threads_map"][thread]["src_local_t"] = np.s_[0:self.parallel_map[self.rank]["threads_map"][thread]["src_local_t_size"]]

            # OVERLAP
            self.parallel_map[self.rank]["threads_map"][thread]["src_local_t_size_overlap"] = self.parallel_map[self.rank]["threads_map"][thread][
                "src_global_t_size_overlap"]

            self.parallel_map[self.rank]["threads_map"][thread]["src_local_t_overlap"] = np.s_[
                0:self.parallel_map[self.rank]["threads_map"][thread]["src_local_t_size_overlap"]]

    # Axis
    def find_time_index(self, t, method="fast", domain="source"):
        """Retourne l'index de la date la plus proche à TIME_DELTA_MIN prêt.
    @type t: datetime ou int
    @param t: date souhaitée ou l'index de la date souhaitée
    @return:  l'index de la date la plus proche à TIME_DELTA_MIN prêt ou une erreur si aucune date n'a pu être trouvée."""

        if type(t) == int or type(t) == np.int32 or type(t) == np.int64:

            if t < 0 or t >= self.get_t_size():
                raise ValueError("Time index have to range between 0 and " + str(
                    self.get_t_size() - 1) + ". Actually Time index = " + str(t))

            return t;

        if type(t) == datetime or type(t) == cftime._cftime.datetime or type(t) == cftime._cftime.real_datetime:

            target_timestamp = (t - TimeCoverage.TIME_DATUM).total_seconds()
            array = np.asarray(self.read_axis_t(type="source", timestamp=1))

            logging.debug("[TimeCoverage][find_time_index()] Looking for : " + str(t))

            if method == "fast":

                nearest_t_index = (np.abs(array - target_timestamp)).argmin()
                if target_timestamp - array[nearest_t_index] == 0.0 or abs(
                        target_timestamp - array[nearest_t_index]) < (TimeCoverage.TIME_DELTA).total_seconds():
                    logging.debug("[TimeCoverage][find_time_index()] Nearest datetime found : " + str(
                        self.read_axis_t(type="source", timestamp=0)[nearest_t_index]))

                    if domain == "source":
                        return nearest_t_index
                    elif domain == "source_global":
                        return self.parallel_map[self.rank]["src_global_t"].start + nearest_t_index
                    else:
                        raise ValueError("Type doesn't match [source, source_global]")
            else:
                raise NotImplementedError("Method " + str(method) + " is not implemented for regular grid.")

            raise NotFoundInRankError(self.rank,
                                      "'" + str(t) + "' not found. Maybe the TimeCoverage.TIME_DELTA (" + str(
                                          TimeCoverage.TIME_DELTA) + ") is too small or the date is out the range.")

        else:
            raise ValueError("" + str(t) + " have to be an integer or a datetime. Current type: " + str(type(t)))


    def __read_thread_axis_t(self,current_thread: int, type="target", with_overlap=False, timestamp=0):
        """Retourne les valeurs de l'axe t.
    @param timestamp: égale 1 si le temps est souhaité en timestamp depuis TIME_DATUM.
    @return:  un tableau à une dimensions [z] au format datetime ou timestamp si timestamp=1."""

        if type == "source" and with_overlap is True:
            return self.reader.read_axis_t(self.parallel_map[self.rank]["threads_map"][current_thread]["src_global_t_overlap"].start,
                                           self.parallel_map[self.rank]["threads_map"][current_thread]["src_global_t_overlap"].stop, timestamp)

        elif type == "source" and with_overlap is False:
            return self.reader.read_axis_t(self.parallel_map[self.rank]["threads_map"][current_thread]["src_global_t"].start,
                                           self.parallel_map[self.rank]["threads_map"][current_thread]["src_global_t"].stop, timestamp)

        elif type == "target" and with_overlap is True:
            if timestamp == 1:
                return [(t - TimeCoverage.TIME_DATUM).total_seconds() \
                        for t in self.target_global_axis_t[self.parallel_map[self.rank]["threads_map"][current_thread]["dst_global_t_overlap"]]];
            return self.target_global_axis_t[self.parallel_map[self.rank]["threads_map"][current_thread]["dst_global_t_overlap"]]

        elif type == "target" and with_overlap is False:
            if timestamp == 1:
                return [(t - TimeCoverage.TIME_DATUM).total_seconds() \
                        for t in self.target_global_axis_t[self.parallel_map[self.rank]["threads_map"][current_thread]["dst_global_t"]]];
            return self.target_global_axis_t[self.parallel_map[self.rank]["threads_map"][current_thread]["dst_global_t"]]

    def read_axis_t(self, type="target", with_overlap=False, timestamp=0):
        """Retourne les valeurs de l'axe t.
    @param timestamp: égale 1 si le temps est souhaité en timestamp depuis TIME_DATUM.
    @return:  un tableau à une dimensions [z] au format datetime ou timestamp si timestamp=1."""
        if type == "target_global":
            if timestamp == 1:
                return [(t - TimeCoverage.TIME_DATUM).total_seconds() \
                        for t in self.target_global_axis_t];
            return self.target_global_axis_t

        elif type == "source_global":
            if timestamp == 1:
                return [(t - TimeCoverage.TIME_DATUM).total_seconds() \
                        for t in self.source_global_axis_t];
            return self.source_global_axis_t

        elif type == "source" and with_overlap is True:
            return self.reader.read_axis_t(self.parallel_map[self.rank]["src_global_t_overlap"].start,
                                           self.parallel_map[self.rank]["src_global_t_overlap"].stop, timestamp)

        elif type == "source" and with_overlap is False:
            return self.reader.read_axis_t(self.parallel_map[self.rank]["src_global_t"].start,
                                           self.parallel_map[self.rank]["src_global_t"].stop, timestamp)

        elif type == "target" and with_overlap is True:
            if timestamp == 1:
                return [(t - TimeCoverage.TIME_DATUM).total_seconds() \
                        for t in self.target_global_axis_t[self.parallel_map[self.rank]["dst_global_t_overlap"]]];
            return self.target_global_axis_t[self.parallel_map[self.rank]["dst_global_t_overlap"]]

        else:
            if timestamp == 1:
                return [(t - TimeCoverage.TIME_DATUM).total_seconds() \
                        for t in self.target_global_axis_t[self.parallel_map[self.rank]["dst_global_t"]]];
            return self.target_global_axis_t[self.parallel_map[self.rank]["dst_global_t"]]

    def get_t_size(self, type="target", with_overlap=False):
        if type == "target_global":
            return self.target_global_t_size
        elif type == "source_global":
            return self.source_global_t_size
        elif type == "source":
            return self.parallel_map[self.rank]["src_local_t_size"]
        elif type == "target" and with_overlap is True:
            return self.parallel_map[self.rank]["dst_local_t_size_overlap"]
        else:
            return self.parallel_map[self.rank]["dst_local_t_size"]

    def __read_variable(self, function_name, time):

        fn = getattr(self.reader, function_name)

        index_t = self.find_time_index(time);

        if self.horizontal_resampling:

            # We use multithreading to compute resampling in parallel
            local_data = np.zeros([self.get_y_size(), self.get_x_size()])
            local_data[:] = np.nan

            with ProcessPoolExecutor(max_workers=self.threads_number) as executor:
                futures = []
                for current_thread in range(0, executor._max_workers):
                    data = fn(
                        self.parallel_map[self.rank]["src_global_t"].start + index_t,
                        self.parallel_map[self.rank]["threads_map"][current_thread]["src_global_x_overlap"].start,
                        self.parallel_map[self.rank]["threads_map"][current_thread]["src_global_x_overlap"].stop,
                        self.parallel_map[self.rank]["threads_map"][current_thread]["src_global_y_overlap"].start,
                        self.parallel_map[self.rank]["threads_map"][current_thread]["src_global_y_overlap"].stop)

                    futures.append(executor.submit(resample_faster_2d_to_grid,
                                                   self.tri[self.rank][current_thread],
                                                   self.read_thread_axis_x(type="target", with_overlap=True,
                                                                           current_thread=current_thread),
                                                   self.read_thread_axis_y(type="target", with_overlap=True,
                                                                           current_thread=current_thread),
                                                   data,
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

            data = fn(
                self.parallel_map[self.rank]["src_global_t"].start + index_t,
                self.parallel_map[self.rank]["src_global_x_overlap"].start,
                self.parallel_map[self.rank]["src_global_x_overlap"].stop,
                self.parallel_map[self.rank]["src_global_y_overlap"].start,
                self.parallel_map[self.rank]["src_global_y_overlap"].stop)

            return data[self.parallel_map[self.rank]["dst_local_y"], self.parallel_map[self.rank]["dst_local_x"]]

    # Variables
    def read_variable_2D_sea_binary_mask_at_time(self, t):
        """Retourne le masque à la date souhaitée sur toute la couverture horizontale.
    @type t: datetime ou l'index
    @param t: date souhaitée
    @return: un tableau en deux dimensions [y,x]."""
        return self.__read_variable(inspect.stack()[0][3],time=t)



    def read_variable_2D_wet_binary_mask_at_time(self, t):
        """Retourne le masque à la date souhaitée sur toute la couverture horizontale.
    @type t: datetime ou l'index
    @param t: date souhaitée
    @return: un tableau en deux dimensions [y,x]."""
        return self.__read_variable(inspect.stack()[0][3],time=t)

    def read_variable_2D_land_binary_mask_at_time(self, t):
        return self.__read_variable(inspect.stack()[0][3],time=t)

    #################
    # HYDRO
    # Sea Surface
    #################
    def read_variable_sea_surface_height_above_mean_sea_level_at_time(self, t):
        return self.__read_variable(inspect.stack()[0][3], time=t)

    def read_variable_sea_surface_height_above_geoid_at_time(self, t):
        return self.__read_variable(inspect.stack()[0][3], time=t)

    def read_variable_sea_water_column_thickness_at_time(self, t):
        return self.__read_variable(inspect.stack()[0][3], time=t)

    def read_variable_sea_surface_temperature_at_time(self, t):
        """Retourne la temperature de surface à la date souhaitée
    @type t: datetime ou l'index
    @param t: date souhaitée
    @return: un tableau en deux dimensions [y,x]."""
        return self.__read_variable(inspect.stack()[0][3], time=t)

    def read_variable_sea_surface_salinity_at_time(self, t):
        """Retourne la salinité de surface à la date souhaitée
    @type t: datetime ou l'index
    @param t: date souhaitée
    @return: un tableau en deux dimensions [y,x]."""
        return self.__read_variable(inspect.stack()[0][3], time=t)

    def read_variable_sea_surface_pressure_at_time(self, t):
        """Retourne la pression à la surface de la mer (sea surface pressure) à la date souhaitée sur toute la couverture horizontale.
    @type t: datetime ou l'index
    @param t: date souhaitée
    @return: un tableau en deux dimensions [y,x]."""
        return self.__read_variable(inspect.stack()[0][3], time=t)

    def read_variable_sea_surface_density_at_time(self, t):
        """Retourne la densité de l'eau de surface à la date souhaitée
    @type t: datetime ou l'index
    @param t: date souhaitée
    @return: un tableau en deux dimensions [y,x]."""

        return self.__read_variable(inspect.stack()[0][3],time=t)

    def read_variable_sea_water_turbidity_at_time(self, t):
        """Retourne la turbidité de l'eau de surface à la date souhaitée
    @type t: datetime ou l'index
    @param t: date souhaitée
    @return: un tableau en deux dimensions [y,x]."""

        return self.__read_variable(inspect.stack()[0][3], time=t)

    def read_variable_sea_water_velocity_at_sea_water_surface_at_time(self, t):
        """Retourne les composantes u,v du courant à la date souhaitée
    @type t: datetime ou l'index
    @param t: date souhaitée
    @return: un tableau en deux dimensions [u_comp,v_comp] contenant chacun deux dimensions [y,x]."""

        index_t = self.find_time_index(t);

        data = self.reader.read_variable_sea_water_velocity_at_sea_water_surface_at_time(
            self.parallel_map[self.rank]["src_global_t"].start + index_t,
            self.parallel_map[self.rank]["src_global_x_overlap"].start,
            self.parallel_map[self.rank]["src_global_x_overlap"].stop,
            self.parallel_map[self.rank]["src_global_y_overlap"].start,
            self.parallel_map[self.rank]["src_global_y_overlap"].stop)

        if self.horizontal_resampling:
            return resample_2d_to_grid(self.read_axis_x(type="source", with_overlap=True),
                                       self.read_axis_y(type="source", with_overlap=True),
                                       self.read_axis_x(type="target", with_overlap=True),
                                       self.read_axis_y(type="target", with_overlap=True),
                                       data[0],
                                       Coverage.HORIZONTAL_INTERPOLATION_METHOD)[
                self.parallel_map[self.rank]["dst_local_y"], self.parallel_map[self.rank]["dst_local_x"]], \
                resample_2d_to_grid(self.read_axis_x(type="source", with_overlap=True),
                                    self.read_axis_y(type="source", with_overlap=True),
                                    self.read_axis_x(type="target", with_overlap=True),
                                    self.read_axis_y(type="target", with_overlap=True),
                                    data[1],
                                    Coverage.HORIZONTAL_INTERPOLATION_METHOD)[
                    self.parallel_map[self.rank]["dst_local_y"], self.parallel_map[self.rank]["dst_local_x"]]

        return data[0][self.parallel_map[self.rank]["dst_local_y"], self.parallel_map[self.rank]["dst_local_x"]], data[1][
            self.parallel_map[self.rank]["dst_local_y"], self.parallel_map[self.rank]["dst_local_x"]]

    #################
    # HYDRO
    # Ground level
    #################

    def read_variable_sea_water_temperature_at_ground_level_at_time(self, t):
        """Retourne la temperature de surface à la date souhaitée
    @type t: datetime ou l'index
    @param t: date souhaitée
    @return: un tableau en deux dimensions [y,x]."""
        return self.__read_variable(inspect.stack()[0][3], time=t)

    def read_variable_sea_water_salinity_at_ground_level_at_time(self, t):
        """Retourne la salinité de surface à la date souhaitée
    @type t: datetime ou l'index
    @param t: date souhaitée
    @return: un tableau en deux dimensions [y,x]."""
        return self.__read_variable(inspect.stack()[0][3], time=t)

    def read_variable_sea_water_velocity_at_ground_level_at_time(self, t):
        """Retourne les composantes u,v du courant à la date souhaitée
           @type t: datetime ou l'index
           @param t: date souhaitée
           @return: un tableau en deux dimensions [u_comp,v_comp] contenant chacun deux dimensions [y,x]."""

        index_t = self.find_time_index(t);

        data = self.reader.read_variable_sea_water_velocity_at_ground_level_at_time(
            self.parallel_map[self.rank]["src_global_t"].start + index_t,
            self.parallel_map[self.rank]["src_global_x_overlap"].start,
            self.parallel_map[self.rank]["src_global_x_overlap"].stop,
            self.parallel_map[self.rank]["src_global_y_overlap"].start,
            self.parallel_map[self.rank]["src_global_y_overlap"].stop)

        if self.horizontal_resampling:
            return resample_2d_to_grid(self.read_axis_x(type="source", with_overlap=True),
                                       self.read_axis_y(type="source", with_overlap=True),
                                       self.read_axis_x(type="target", with_overlap=True),
                                       self.read_axis_y(type="target", with_overlap=True),
                                       data[0],
                                       Coverage.HORIZONTAL_INTERPOLATION_METHOD)[
                self.parallel_map[self.rank]["dst_local_y"], self.parallel_map[self.rank]["dst_local_x"]], \
                resample_2d_to_grid(self.read_axis_x(type="source", with_overlap=True),
                                    self.read_axis_y(type="source", with_overlap=True),
                                    self.read_axis_x(type="target", with_overlap=True),
                                    self.read_axis_y(type="target", with_overlap=True),
                                    data[1],
                                    Coverage.HORIZONTAL_INTERPOLATION_METHOD)[
                    self.parallel_map[self.rank]["dst_local_y"], self.parallel_map[self.rank]["dst_local_x"]]

        return data[0][self.parallel_map[self.rank]["dst_local_y"], self.parallel_map[self.rank]["dst_local_x"]], data[1][
            self.parallel_map[self.rank]["dst_local_y"], self.parallel_map[self.rank]["dst_local_x"]]

    #################
    # HYDRO
    # 2D
    #################
    def read_variable_barotropic_sea_water_velocity_at_time(self, t):
        """Retourne les composantes u,v du courant à la date souhaitée
    @type t: datetime ou l'index
    @param t: date souhaitée
    @return: un tableau en deux dimensions [u_comp,v_comp] contenant chacun deux dimensions [y,x]."""

        index_t = self.find_time_index(t);

        data = self.reader.read_variable_barotropic_sea_water_velocity_at_time(
            self.parallel_map[self.rank]["src_global_t"].start + index_t,
            self.parallel_map[self.rank]["src_global_x_overlap"].start,
            self.parallel_map[self.rank]["src_global_x_overlap"].stop,
            self.parallel_map[self.rank]["src_global_y_overlap"].start,
            self.parallel_map[self.rank]["src_global_y_overlap"].stop)

        if self.horizontal_resampling:
            return resample_2d_to_grid(self.read_axis_x(type="source", with_overlap=True),
                                       self.read_axis_y(type="source", with_overlap=True),
                                       self.read_axis_x(type="target", with_overlap=True),
                                       self.read_axis_y(type="target", with_overlap=True),
                                       data[0],
                                       Coverage.HORIZONTAL_INTERPOLATION_METHOD)[
                self.parallel_map[self.rank]["dst_local_y"], self.parallel_map[self.rank]["dst_local_x"]], \
                resample_2d_to_grid(self.read_axis_x(type="source", with_overlap=True),
                                    self.read_axis_y(type="source", with_overlap=True),
                                    self.read_axis_x(type="target", with_overlap=True),
                                    self.read_axis_y(type="target", with_overlap=True),
                                    data[1],
                                    Coverage.HORIZONTAL_INTERPOLATION_METHOD)[
                    self.parallel_map[self.rank]["dst_local_y"], self.parallel_map[self.rank]["dst_local_x"]]

        return data[0][self.parallel_map[self.rank]["dst_local_y"], self.parallel_map[self.rank]["dst_local_x"]], data[1][
            self.parallel_map[self.rank]["dst_local_y"], self.parallel_map[self.rank]["dst_local_x"]]

    def read_variable_barotropic_sea_water_speed_at_time(self, date):
        comp = self.read_variable_barotropic_sea_water_velocity_at_time(date)
        result = np.zeros([self.get_y_size(), self.get_x_size()])
        result[:] = np.nan
        for x in range(0, self.get_x_size()):
            for y in range(0, self.get_y_size()):
                result[y, x] = math.sqrt(comp[0][y, x] ** 2 + comp[1][y, x] ** 2)

        return result

    def read_variable_barotropic_sea_water_from_direction_at_time(self, date):
        comp = self.read_variable_barotropic_sea_water_velocity_at_time(date)
        result = np.zeros([self.get_y_size(), self.get_x_size()])
        result[:] = np.nan
        for x in range(0, self.get_x_size()):
            for y in range(0, self.get_y_size()):
                result[y, x] = 270. - (180.0 / math.pi) * (math.atan2(comp[0][y, x], comp[1][y, x])) + 180.0 % 360.0

        return result

    def read_variable_barotropic_sea_water_to_direction_at_time(self, date):
        comp = self.read_variable_barotropic_sea_water_velocity_at_time(date)
        result = np.zeros([self.get_y_size(), self.get_x_size()])
        result[:] = np.nan
        for x in range(0, self.get_x_size()):
            for y in range(0, self.get_y_size()):
                result[y, x] = 270. - (180.0 / math.pi) * (math.atan2(comp[0][y, x], comp[1][y, x])) % 360.0

        return result

    #################
    # WAVES
    # Sea Surface
    #################
    def read_variable_sea_surface_wave_significant_height_at_time(self, t):
        """Retourne la hauteur significative des vagues à la date souhaitée sur toute la couverture horizontale.
    @type t: datetime ou l'index
    @param t: date souhaitée
    @return: un tableau en deux dimensions [y,x]."""

        return self.__read_variable(inspect.stack()[0][3], time=t)

    def read_variable_sea_surface_wave_breaking_height_at_time(self, t):
        """Retourne la hauteur de déferlement des vagues à la date souhaitée sur toute la couverture horizontale.
    @type t: datetime ou l'index
    @param t: date souhaitée
    @return: un tableau en deux dimensions [y,x]."""
        return self.__read_variable(inspect.stack()[0][3], time=t)

    def read_variable_sea_surface_wave_mean_period_at_time(self, t):
        return self.__read_variable(inspect.stack()[0][3], time=t)

    def read_variable_sea_surface_wave_peak_period_at_time(self, t):

        return self.__read_variable(inspect.stack()[0][3], time=t)

    def read_variable_sea_surface_wave_from_direction_at_time(self, t):

        return self.__read_variable(inspect.stack()[0][3],time=t)

    def read_variable_sea_surface_wave_to_direction_at_time(self, t):

        return self.__read_variable(inspect.stack()[0][3],time=t)

    def read_variable_sea_surface_wave_stokes_drift_velocity_at_time(self, t):
        """Retourne la dérive de Stokes en surface à la date souhaitée sur toute la couverture horizontale.
    @type t: datetime ou l'index
    @param t: date souhaitée
    @return: un tableau en deux dimensions [y,x]."""

        index_t = self.find_time_index(t);

        data = self.reader.read_variable_sea_surface_wave_stokes_drift_velocity_at_time(
            self.parallel_map[self.rank]["src_global_t"].start + index_t,
            self.parallel_map[self.rank]["src_global_x_overlap"].start,
            self.parallel_map[self.rank]["src_global_x_overlap"].stop,
            self.parallel_map[self.rank]["src_global_y_overlap"].start,
            self.parallel_map[self.rank]["src_global_y_overlap"].stop)

        if self.horizontal_resampling:
            return resample_2d_to_grid(self.read_axis_x(type="source", with_overlap=True),
                                       self.read_axis_y(type="source", with_overlap=True),
                                       self.read_axis_x(type="target", with_overlap=True),
                                       self.read_axis_y(type="target", with_overlap=True),
                                       data[0],
                                       Coverage.HORIZONTAL_INTERPOLATION_METHOD)[
                self.parallel_map[self.rank]["dst_local_y"], self.parallel_map[self.rank]["dst_local_x"]], \
                resample_2d_to_grid(self.read_axis_x(type="source", with_overlap=True),
                                    self.read_axis_y(type="source", with_overlap=True),
                                    self.read_axis_x(type="target", with_overlap=True),
                                    self.read_axis_y(type="target", with_overlap=True),
                                    data[1],
                                    Coverage.HORIZONTAL_INTERPOLATION_METHOD)[
                    self.parallel_map[self.rank]["dst_local_y"], self.parallel_map[self.rank]["dst_local_x"]]

        return data[0][self.parallel_map[self.rank]["dst_local_y"], self.parallel_map[self.rank]["dst_local_x"]], data[1][
            self.parallel_map[self.rank]["dst_local_y"], self.parallel_map[self.rank]["dst_local_x"]]

    def read_variable_radiation_pressure_bernouilli_head_at_time(self, t):
        """Retourne la pression J due aux vagues à la date souhaitée sur toute la couverture horizontale.
    @type t: datetime ou l'index
    @param t: date souhaitée
    @return: un tableau en deux dimensions [y,x]."""

        return self.__read_variable(inspect.stack()[0][3], time=t)

    def read_variable_sea_surface_wave_energy_flux_to_ocean_at_time(self, t):
        """Retourne la waves_to_ocean_energy_flux à la date souhaitée sur toute la couverture horizontale.
    @type t: datetime ou l'index
    @param t: date souhaitée
    @return: un tableau en deux dimensions [y,x]."""

        index_t = self.find_time_index(t);

        data = self.reader.read_variable_sea_surface_wave_energy_flux_to_ocean_at_time(
            self.parallel_map[self.rank]["src_global_t"].start + index_t,
            self.parallel_map[self.rank]["src_global_x_overlap"].start,
            self.parallel_map[self.rank]["src_global_x_overlap"].stop,
            self.parallel_map[self.rank]["src_global_y_overlap"].start,
            self.parallel_map[self.rank]["src_global_y_overlap"].stop)

        if self.horizontal_resampling:
            return resample_2d_to_grid(self.read_axis_x(type="source", with_overlap=True),
                                       self.read_axis_y(type="source", with_overlap=True),
                                       self.read_axis_x(type="target", with_overlap=True),
                                       self.read_axis_y(type="target", with_overlap=True),
                                       data[0],
                                       Coverage.HORIZONTAL_INTERPOLATION_METHOD)[
                self.parallel_map[self.rank]["dst_local_y"], self.parallel_map[self.rank]["dst_local_x"]], \
                resample_2d_to_grid(self.read_axis_x(type="source", with_overlap=True),
                                    self.read_axis_y(type="source", with_overlap=True),
                                    self.read_axis_x(type="target", with_overlap=True),
                                    self.read_axis_y(type="target", with_overlap=True),
                                    data[1],
                                    Coverage.HORIZONTAL_INTERPOLATION_METHOD)[
                    self.parallel_map[self.rank]["dst_local_y"], self.parallel_map[self.rank]["dst_local_x"]]

        return data[0][self.parallel_map[self.rank]["dst_local_y"], self.parallel_map[self.rank]["dst_local_x"]], data[1][
            self.parallel_map[self.rank]["dst_local_y"], self.parallel_map[self.rank]["dst_local_x"]]

    def read_variable_sea_surface_wave_energy_dissipation_at_ground_level_at_time(self, t):
        """Retourne la l'énergie des vagues dissipée par le fond à la date souhaitée sur toute la couverture horizontale.
    @type t: datetime ou l'index
    @param t: date souhaitée
    @return: un tableau en deux dimensions [y,x]."""

        return self.__read_variable(inspect.stack()[0][3],time=t)

    #################
    # WAVES
    # Momentum flux
    #################
    def read_variable_atmosphere_momentum_flux_to_waves_at_time(self, t):
        """Retourne la composante u du tau atmosphere->vagues à la date souhaitée sur toute la couverture horizontale.
    @type t: datetime ou l'index
    @param t: date souhaitée
    @return: un tableau en deux dimensions [y,x]."""

        index_t = self.find_time_index(t);

        data = self.reader.read_variable_atmosphere_momentum_flux_to_waves_at_time(
            self.parallel_map[self.rank]["src_global_t"].start + index_t,
            self.parallel_map[self.rank]["src_global_x_overlap"].start,
            self.parallel_map[self.rank]["src_global_x_overlap"].stop,
            self.parallel_map[self.rank]["src_global_y_overlap"].start,
            self.parallel_map[self.rank]["src_global_y_overlap"].stop)

        if self.horizontal_resampling:
            return resample_2d_to_grid(self.read_axis_x(type="source", with_overlap=True),
                                       self.read_axis_y(type="source", with_overlap=True),
                                       self.read_axis_x(type="target", with_overlap=True),
                                       self.read_axis_y(type="target", with_overlap=True),
                                       data[0],
                                       Coverage.HORIZONTAL_INTERPOLATION_METHOD)[
                self.parallel_map[self.rank]["dst_local_y"], self.parallel_map[self.rank]["dst_local_x"]], \
                resample_2d_to_grid(self.read_axis_x(type="source", with_overlap=True),
                                    self.read_axis_y(type="source", with_overlap=True),
                                    self.read_axis_x(type="target", with_overlap=True),
                                    self.read_axis_y(type="target", with_overlap=True),
                                    data[1],
                                    Coverage.HORIZONTAL_INTERPOLATION_METHOD)[
                    self.parallel_map[self.rank]["dst_local_y"], self.parallel_map[self.rank]["dst_local_x"]]

        return data[0][self.parallel_map[self.rank]["dst_local_y"], self.parallel_map[self.rank]["dst_local_x"]], data[1][
            self.parallel_map[self.rank]["dst_local_y"], self.parallel_map[self.rank]["dst_local_x"]]

    def read_variable_waves_momentum_flux_to_ocean_at_time(self, t):
        """Retourne la composante u du tau vagues->ocean à la date souhaitée sur toute la couverture horizontale.
    @type t: datetime ou l'index
    @param t: date souhaitée
    @return: un tableau en deux dimensions [y,x]."""

        index_t = self.find_time_index(t);

        data = self.reader.read_variable_waves_momentum_flux_to_ocean_at_time(
            self.parallel_map[self.rank]["src_global_t"].start + index_t,
            self.parallel_map[self.rank]["src_global_x_overlap"].start,
            self.parallel_map[self.rank]["src_global_x_overlap"].stop,
            self.parallel_map[self.rank]["src_global_y_overlap"].start,
            self.parallel_map[self.rank]["src_global_y_overlap"].stop)

        if self.horizontal_resampling:
            return resample_2d_to_grid(self.read_axis_x(type="source", with_overlap=True),
                                       self.read_axis_y(type="source", with_overlap=True),
                                       self.read_axis_x(type="target", with_overlap=True),
                                       self.read_axis_y(type="target", with_overlap=True),
                                       data[0],
                                       Coverage.HORIZONTAL_INTERPOLATION_METHOD)[
                self.parallel_map[self.rank]["dst_local_y"], self.parallel_map[self.rank]["dst_local_x"]], \
                resample_2d_to_grid(self.read_axis_x(type="source", with_overlap=True),
                                    self.read_axis_y(type="source", with_overlap=True),
                                    self.read_axis_x(type="target", with_overlap=True),
                                    self.read_axis_y(type="target", with_overlap=True),
                                    data[1],
                                    Coverage.HORIZONTAL_INTERPOLATION_METHOD)[
                    self.parallel_map[self.rank]["dst_local_y"], self.parallel_map[self.rank]["dst_local_x"]]

        return data[0][self.parallel_map[self.rank]["dst_local_y"], self.parallel_map[self.rank]["dst_local_x"]], data[1][
            self.parallel_map[self.rank]["dst_local_y"], self.parallel_map[self.rank]["dst_local_x"]]

    #################
    # METEO
    # 2D
    #################
    def read_variable_rainfall_amount_at_time(self, t):
        """Retourne les composantes u,v de rain à la date souhaitée
    @type t: datetime ou l'index
    @param t: date souhaitée
    @return: un tableau en deux dimensions [u_comp,v_comp] contenant chacun deux dimensions [y,x]."""

        return self.__read_variable(inspect.stack()[0][3],time=t)

    #################
    # METEO
    # Sea surface
    #################

    def read_variable_surface_air_pressure_at_time(self, t):
        """Retourne la pression à la surface à la date souhaitée sur toute la couverture horizontale.
    @type t: datetime ou l'index
    @param t: date souhaitée
    @return: un tableau en deux dimensions [y,x]."""

        return self.__read_variable(inspect.stack()[0][3],time=t)

    def read_variable_sea_surface_air_pressure_at_time(self, t):
        """Retourne la pression à la surface à la date souhaitée sur toute la couverture horizontale.
    @type t: datetime ou l'index
    @param t: date souhaitée
    @return: un tableau en deux dimensions [y,x]."""

        return self.__read_variable(inspect.stack()[0][3], time=t)

    def read_variable_wind_stress_at_time(self, t):
        """Retourne les composantes u,v de la contrainte de vent à la date souhaitée
    @type t: datetime ou l'index
    @param t: date souhaitée
    @return: un tableau en deux dimensions [u_comp,v_comp] contenant chacun deux dimensions [y,x]."""

        index_t = self.find_time_index(t);

        data = self.reader.read_variable_wind_stress_at_time(
            self.parallel_map[self.rank]["src_global_t"].start + index_t,
            self.parallel_map[self.rank]["src_global_x_overlap"].start,
            self.parallel_map[self.rank]["src_global_x_overlap"].stop,
            self.parallel_map[self.rank]["src_global_y_overlap"].start,
            self.parallel_map[self.rank]["src_global_y_overlap"].stop)

        if self.horizontal_resampling:
            return resample_2d_to_grid(self.read_axis_x(type="source", with_overlap=True),
                                       self.read_axis_y(type="source", with_overlap=True),
                                       self.read_axis_x(type="target", with_overlap=True),
                                       self.read_axis_y(type="target", with_overlap=True),
                                       data[0],
                                       Coverage.HORIZONTAL_INTERPOLATION_METHOD)[
                self.parallel_map[self.rank]["dst_local_y"], self.parallel_map[self.rank]["dst_local_x"]], \
                resample_2d_to_grid(self.read_axis_x(type="source", with_overlap=True),
                                    self.read_axis_y(type="source", with_overlap=True),
                                    self.read_axis_x(type="target", with_overlap=True),
                                    self.read_axis_y(type="target", with_overlap=True),
                                    data[1],
                                    Coverage.HORIZONTAL_INTERPOLATION_METHOD)[
                    self.parallel_map[self.rank]["dst_local_y"], self.parallel_map[self.rank]["dst_local_x"]]

        return data[0][self.parallel_map[self.rank]["dst_local_y"], self.parallel_map[self.rank]["dst_local_x"]], data[1][
            self.parallel_map[self.rank]["dst_local_y"], self.parallel_map[self.rank]["dst_local_x"]]

    def read_variable_surface_downward_sensible_heat_flux_at_time(self, t):
        """Retourne les composantes u,v de surface sensible heat flux à la date souhaitée
    @type t: datetime ou l'index
    @param t: date souhaitée
    @return: un tableau en deux dimensions [u_comp,v_comp] contenant chacun deux dimensions [y,x]."""

        return self.__read_variable(inspect.stack()[0][3], time=t)

    def read_variable_surface_downward_latent_heat_flux_at_time(self, t):
        """Retourne les composantes u,v de surface latente heat flux à la date souhaitée
    @type t: datetime ou l'index
    @param t: date souhaitée
    @return: un tableau en deux dimensions [u_comp,v_comp] contenant chacun deux dimensions [y,x]."""

        return self.__read_variable(inspect.stack()[0][3], time=t)

    def read_variable_surface_air_temperature_at_time(self, t):
        """Retourne les composantes u,v de surface air temperature à la date souhaitée
    @type t: datetime ou l'index
    @param t: date souhaitée
    @return: un tableau en deux dimensions [u_comp,v_comp] contenant chacun deux dimensions [y,x]."""

        return self.__read_variable(inspect.stack()[0][3], time=t)

    def read_variable_dew_point_temperature_at_time(self, t):
        """Retourne les composantes u,v de dewpoint temperature à la date souhaitée
    @type t: datetime ou l'index
    @param t: date souhaitée
    @return: un tableau en deux dimensions [u_comp,v_comp] contenant chacun deux dimensions [y,x]."""

        return self.__read_variable(inspect.stack()[0][3], time=t)

    def read_variable_surface_downward_solar_radiation_at_time(self, t):
        """Retourne les composantes u,v de surface solar radiation downwards à la date souhaitée
    @type t: datetime ou l'index
    @param t: date souhaitée
    @return: un tableau en deux dimensions [u_comp,v_comp] contenant chacun deux dimensions [y,x]."""

        return self.__read_variable(inspect.stack()[0][3],time=t)

    def read_variable_surface_downward_thermal_radiation_at_time(self, t):
        """Retourne les composantes u,v de surface thermal radiation downwards à la date souhaitée
    @type t: datetime ou l'index
    @param t: date souhaitée
    @return: un tableau en deux dimensions [u_comp,v_comp] contenant chacun deux dimensions [y,x]."""

        return self.__read_variable(inspect.stack()[0][3], time=t)

    def read_variable_surface_solar_radiation_at_time(self, t):
        """Retourne les composantes u,v de surface solar radiation à la date souhaitée
    @type t: datetime ou l'index
    @param t: date souhaitée
    @return: un tableau en deux dimensions [u_comp,v_comp] contenant chacun deux dimensions [y,x]."""

        return self.__read_variable(inspect.stack()[0][3],time=t)

    def read_variable_surface_thermal_radiation_at_time(self, t):
        """Retourne les composantes u,v de surface thermal radiation à la date souhaitée
    @type t: datetime ou l'index
    @param t: date souhaitée
    @return: un tableau en deux dimensions [u_comp,v_comp] contenant chacun deux dimensions [y,x]."""

        return self.__read_variable(inspect.stack()[0][3], time=t)

    #################
    # METEO
    # At 10 m
    #################
    def read_variable_wind_10m_at_time(self, t):
        """Retourne les composantes u,v du vent à la date souhaitée
    @type t: datetime ou l'index
    @param t: date souhaitée
    @return: un tableau en deux dimensions [u_comp,v_comp] contenant chacun deux dimensions [y,x]."""

        index_t = self.find_time_index(t);

        data = self.reader.read_variable_wind_10m_at_time(
            self.parallel_map[self.rank]["src_global_t"].start + index_t,
            self.parallel_map[self.rank]["src_global_x_overlap"].start,
            self.parallel_map[self.rank]["src_global_x_overlap"].stop,
            self.parallel_map[self.rank]["src_global_y_overlap"].start,
            self.parallel_map[self.rank]["src_global_y_overlap"].stop)

        if self.horizontal_resampling:
            return resample_2d_to_grid(self.read_axis_x(type="source", with_overlap=True),
                                       self.read_axis_y(type="source", with_overlap=True),
                                       self.read_axis_x(type="target", with_overlap=True),
                                       self.read_axis_y(type="target", with_overlap=True),
                                       data[0],
                                       Coverage.HORIZONTAL_INTERPOLATION_METHOD)[
                self.parallel_map[self.rank]["dst_local_y"], self.parallel_map[self.rank]["dst_local_x"]], \
                resample_2d_to_grid(self.read_axis_x(type="source", with_overlap=True),
                                    self.read_axis_y(type="source", with_overlap=True),
                                    self.read_axis_x(type="target", with_overlap=True),
                                    self.read_axis_y(type="target", with_overlap=True),
                                    data[1],
                                    Coverage.HORIZONTAL_INTERPOLATION_METHOD)[
                    self.parallel_map[self.rank]["dst_local_y"], self.parallel_map[self.rank]["dst_local_x"]]

        return data[0][self.parallel_map[self.rank]["dst_local_y"], self.parallel_map[self.rank]["dst_local_x"]], data[1][
            self.parallel_map[self.rank]["dst_local_y"], self.parallel_map[self.rank]["dst_local_x"]]

    def read_variable_wind_speed_10m_at_time(self, date):
        comp = self.read_variable_wind_10m_at_time(date)
        result = np.zeros([self.get_y_size(), self.get_x_size()])
        result[:] = np.nan
        for x in range(0, self.get_x_size()):
            for y in range(0, self.get_y_size()):
                result[y, x] = math.sqrt(comp[0][y, x] ** 2 + comp[1][y, x] ** 2)

        return result

    def read_variable_wind_from_direction_10m_at_time(self, date):
        comp = self.read_variable_wind_10m_at_time(date)
        result = np.zeros([self.get_y_size(), self.get_x_size()])
        result[:] = np.nan
        for x in range(0, self.get_x_size()):
            for y in range(0, self.get_y_size()):
                result[y, x] = 270. - (180.0 / math.pi) * (math.atan2(comp[0][y, x], comp[1][y, x])) + 180.0 % 360.0

        return result

    def read_variable_wind_to_direction_10m_at_time(self, date):
        comp = self.read_variable_wind_10m_at_time(date)
        result = np.zeros([self.get_y_size(), self.get_x_size()])
        result[:] = np.nan
        for x in range(0, self.get_x_size()):
            for y in range(0, self.get_y_size()):
                result[y, x] = 270. - (180.0 / math.pi) * (math.atan2(comp[0][y, x], comp[1][y, x])) % 360.0

        return result
