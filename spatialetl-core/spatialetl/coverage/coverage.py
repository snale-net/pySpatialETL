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

import os
from importlib.util import find_spec

import numpy as np
from array_split import shape_split

from spatialetl.exception.not_found_in_rank_error import NotFoundInRankError
from spatialetl.operator.interpolator.interpolator_core import resample_2d_to_grid
from spatialetl.utils.distance import distance_on_unit_sphere
from spatialetl.utils.logger import logging

mpi_lib = find_spec("mpi4py")
MPI_FOUND = mpi_lib is not None


class Coverage(object):
    """
    The Coverage class represents a spatial coverage on the horizontal plane.

    The points representing this coverage can be aligned on a regular grid (x,y) or on an irregular grid ((x1,y1),(x2,y2)).
    Depending on the grid type, the axis reading functions will return one or two-dimensional arrays.
    To avoid loading the entire file into memory, the coverage contains a pointer to a reader.
    Layers are therefore read on demand from the file.

    Note: The axes are always reversed in the arrays because of NetCDF.
    So the y-axis first, then the x-axis. Example: [y,x]

    Examples
    --------
    >>> coverage = Coverage(myReader, bbox=[-5, 5, -5, 5], resolution_x=0.1, resolution_y=0.1)
    """

    HORIZONTAL_INTERPOLATION_METHOD = "linear"
    HORIZONTAL_OVERLAPING_SIZE = 2

    def __init__(self, reader, bbox=None, resolution_x=None, resolution_y=None,parallel=True,nb_thread:int=os.cpu_count()):
        """
        Initialize the Coverage object with a file reader and optional bounding box and resolution.
        
        Parameters
        ----------
        
        reader : object
            File reader instance.
        bbox : list, optional
            Bounding box coordinates [xmin, xmax, ymin, ymax].
        resolution_x : float, optional
            Resolution along the x-axis.
        resolution_y : float, optional
            Resolution along the y-axis.

        Examples
        --------
        >>> coverage = Coverage(reader, bbox=[0, 10, 0, 10], resolution_x=1.0, resolution_y=1.0)
        """

        self.reader = reader;

        if parallel:
            if MPI_FOUND:
                # Parallel MPI
                self.comm = mpi_lib.MPI.COMM_WORLD
                self.size = self.comm.Get_size()
                self.rank = self.comm.Get_rank()
            else:
                # No parallel MPI
                self.map_mpi = None
                self.size = 1
                self.rank = 0

            # Parallel multithreading
            self.threads_number = nb_thread
            self.threading_map = None
        else:
            # No parallel MPI
            self.map_mpi = None
            self.size = 1
            self.rank = 0

            # No parallel multithreading
            self.threads_number = 1
            self.threading_map = None

        self.source_regular_grid = self.reader.is_regular_grid()
        self.target_regular_grid = self.source_regular_grid
        self.horizontal_resampling = False

        self.source_global_x_size = self.reader.get_x_size()
        self.source_global_y_size = self.reader.get_y_size()
        self.source_global_axis_x = self.reader.read_axis_x(0, self.source_global_x_size, 0,
                                                            self.source_global_y_size)
        self.source_global_axis_y = self.reader.read_axis_y(0, self.source_global_x_size, 0,
                                                            self.source_global_y_size)

        # We adjust the computation based on the bbox
        if bbox is None:
            # we compute the destination grid
            Ymin = np.min(self.source_global_axis_y)
            Ymax = np.max(self.source_global_axis_y)
            Xmin = np.min(self.source_global_axis_x)
            Xmax = np.max(self.source_global_axis_x)
        else:
            if self.check_bbox_validity(bbox) is False:
                raise ValueError("Your Bbox is not valid or is outside the coverage")
            Ymin = bbox[2]
            Ymax = bbox[3]
            Xmin = bbox[0]
            Xmax = bbox[1]

        if self.is_regular_grid(type="source"):

            idx = np.where((self.source_global_axis_x >= Xmin) &
                           (self.source_global_axis_x <= Xmax))

            xmin = np.min(idx[0])
            xmax = np.max(idx[0]) + 1

            idx = np.where((self.source_global_axis_y >= Ymin) &
                           (self.source_global_axis_y <= Ymax))

            ymin = np.min(idx[0])
            ymax = np.max(idx[0]) + 1

            self.target_global_axis_x = self.source_global_axis_x[xmin:xmax]
            self.target_global_x_size = xmax - xmin
            self.target_global_axis_y = self.source_global_axis_y[ymin:ymax]
            self.target_global_y_size = ymax - ymin

        else:

            idx = np.where((self.source_global_axis_x >= Xmin) &
                           (self.source_global_axis_x <= Xmax) &
                           (self.source_global_axis_y >= Ymin) &
                           (self.source_global_axis_y <= Ymax))

            if len(idx[0]) == 0:
                raise ValueError("No values found")

            ymin = np.min(idx[0])
            ymax = np.max(idx[0]) + 1
            xmin = np.min(idx[1])
            xmax = np.max(idx[1]) + 1

            self.target_global_axis_x = self.source_global_axis_x[ymin:ymax, xmin:xmax]
            self.target_global_x_size = xmax - xmin
            self.target_global_axis_y = self.source_global_axis_y[ymin:ymax, xmin:xmax]
            self.target_global_y_size = ymax - ymin

        # Compute the destination grid
        self.target_global_res_x = None
        self.target_global_res_y = None

        if resolution_x is not None and resolution_y is not None:

            self.horizontal_resampling = True
            self.target_regular_grid = True

            res = np.mean([resolution_x, resolution_y])
            self.target_global_res_x = res
            self.target_global_res_y = res

            self.target_global_axis_x = np.arange(Xmin, Xmax, res)
            self.target_global_axis_y = np.arange(Ymin, Ymax, res)

            self.target_global_x_size = len(self.target_global_axis_x)
            self.target_global_y_size = len(self.target_global_axis_y)

            if type(self) == Coverage:
                if self.rank == 0:
                    logging.info(
                        '[horizontal_interpolation] Source grid size : (' + str(self.source_global_x_size) + ", " + str(
                            self.source_global_y_size) + ")")
                    logging.info(
                        '[horizontal_interpolation] Target grid size : (' + str(self.target_global_x_size) + ", " + str(
                            self.target_global_y_size) + ")")

        if type(self) == Coverage:
            self.create_mpi_map()
            self.update_mpi_map()

            self.create_threading_map()
            self.update_threading_map()

            if self.rank == 0 and MPI_FOUND:
                logging.debug("MPI map:")
            else:
                logging.debug("Multithreading map:")
                logging.debug(f"{"-" * 10} Target grid {"-" * 10}")
                for key in ['dst_global_x','dst_global_y','dst_local_x_size','dst_local_y_size']:
                    logging.debug(f"{key} = {self.map_mpi[self.rank][key]}")

            if MPI_FOUND:
                logging.debug(f"{"-" * 10} Proc n° {self.rank} {"-" * 10}")
                for key in self.map_mpi[self.rank]:
                    logging.debug(f"    {key} = {self.map_mpi[self.rank][key]}")

            logging.debug(self.threads_number)

            for thread in range(0, self.threads_number):
                logging.debug(f"   {"-" * 10} Thread n° {thread} {"-" * 10}")
                for key in self.map_mpi[self.rank]:
                    logging.debug(f"    {key} = {self.threading_map[self.rank][thread][key]}")

            if self.rank == 0:
                logging.debug("-"* 20)

        # try to fill metadata
        self.read_metadata()

    def check_bbox_validity(self, candidate):
        """
        Check the validity of the bounding box.

        Parameters
        ----------
        candidate : list
            Bounding box coordinates [xmin, xmax, ymin, ymax].

        Returns
        -------
        bool
            True if the bounding box is valid, False otherwise.

        Examples
        --------
        >>> coverage = Coverage(myReader, bbox=[0, 10, 0, 10], resolution_x=1.0, resolution_y=1.0)
        """
        Ymin = candidate[2]
        Ymax = candidate[3]
        Xmin = candidate[0]
        Xmax = candidate[1]

        if Xmax <= Xmin:
            return False
        if Ymax <= Ymin:
            return False

        if Ymin < np.min(self.source_global_axis_y):
            return False
        if Ymax > np.max(self.source_global_axis_y):
            return False
        if Xmin < np.min(self.source_global_axis_x):
            return False
        if Xmax > np.max(self.source_global_axis_x):
            return False

        return True

    def check_point_is_inside(self, target_lon, target_lat, lon, lat, tolerance=5):
        """
        Check if a point is inside the coverage.

        Parameters
        ----------
        target_lon : float
            Longitude of the target point.
        target_lat : float
            Latitude of the target point.
        lon : array
            Longitude array.
        lat : array
            Latitude array.
        tolerance : int, optional
            Tolerance for rounding.

        Returns
        -------
        bool
            True if the point is inside the coverage, False otherwise.

        Examples
        --------
        >>> coverage.check_point_is_inside(5.0, 45.0, lon_array, lat_array)
        True
        """

        if np.round(target_lat, decimals=tolerance) < np.min(np.round(lat, decimals=tolerance)):
            return False
        if np.round(target_lat, decimals=tolerance) > np.max(np.round(lat, decimals=tolerance)):
            return False
        if np.round(target_lon, decimals=tolerance) < np.min(np.round(lon, decimals=tolerance)):
            return False
        if np.round(target_lon, decimals=tolerance) > np.max(np.round(lon, decimals=tolerance)):
            return False

        return True

    def create_mpi_map(self):
        """
        Create the MPI map for parallel processing.
        The MPI map is a dictionary that contains the mapping of the source and destination grids for each MPI rank.

        Examples
        --------
        >>> coverage.create_mpi_map()
        """
        self.map_mpi = np.empty(self.size, dtype=object)
        target_sample = (self.target_global_y_size, self.target_global_x_size)

        # Split the axes
        target_slices = shape_split(target_sample, self.size, axis=[0, 0])

        slice_index = 0
        for slyce in target_slices.flatten():
            slice = tuple(slyce)

            map = {}

            ### Destination grid ###
            map["dst_global_x"] = slice[1]
            map["dst_global_y"] = slice[0]

            map["dst_local_x_size"] = map["dst_global_x"].stop - map["dst_global_x"].start
            map["dst_local_y_size"] = map["dst_global_y"].stop - map["dst_global_y"].start

            dst_global_x_min_overlap = max(0, map["dst_global_x"].start - Coverage.HORIZONTAL_OVERLAPING_SIZE)
            dst_global_x_max_overlap = min(self.target_global_x_size,
                                           map["dst_global_x"].stop + Coverage.HORIZONTAL_OVERLAPING_SIZE)
            map["dst_global_x_overlap"] = np.s_[dst_global_x_min_overlap:dst_global_x_max_overlap]

            dst_global_y_min_overlap = max(0, map["dst_global_y"].start - Coverage.HORIZONTAL_OVERLAPING_SIZE)
            dst_global_y_max_overlap = min(self.target_global_y_size,
                                           map["dst_global_y"].stop + Coverage.HORIZONTAL_OVERLAPING_SIZE)
            map["dst_global_y_overlap"] = np.s_[dst_global_y_min_overlap:dst_global_y_max_overlap]

            map["dst_global_x_size_overlap"] = map["dst_global_x_overlap"].stop - map["dst_global_x_overlap"].start
            map["dst_global_y_size_overlap"] = map["dst_global_y_overlap"].stop - map["dst_global_y_overlap"].start

            dst_x_min = Coverage.HORIZONTAL_OVERLAPING_SIZE
            dst_x_max = map["dst_global_x_size_overlap"] - Coverage.HORIZONTAL_OVERLAPING_SIZE
            dst_y_min = Coverage.HORIZONTAL_OVERLAPING_SIZE
            dst_y_max = map["dst_global_y_size_overlap"] - Coverage.HORIZONTAL_OVERLAPING_SIZE

            if map["dst_global_x"].start == 0:
                dst_x_min = 0

            if map["dst_global_x"].stop == self.target_global_x_size:
                dst_x_max = map["dst_global_x_size_overlap"]

            if map["dst_global_y"].start == 0:
                dst_y_min = 0

            if map["dst_global_y"].stop == self.target_global_y_size:
                dst_y_max = map["dst_global_y_size_overlap"]

            map["dst_local_x"] = np.s_[dst_x_min:dst_x_max]
            map["dst_local_y"] = np.s_[dst_y_min:dst_y_max]

            # Source grille
            map["src_global_x"] = map["dst_global_x"]
            map["src_global_y"] = map["dst_global_y"]

            map["src_global_x_overlap"] = map["dst_global_x_overlap"]
            map["src_global_y_overlap"] = map["dst_global_y_overlap"]

            map["src_local_x"] = map["dst_local_x"]
            map["src_local_y"] = map["dst_local_y"]

            map["src_local_x_size"] = map["dst_local_x_size"]
            map["src_local_y_size"] = map["dst_local_y_size"]

            map["src_local_x_size_overlap"] = map["dst_global_x_size_overlap"]
            map["src_local_y_size_overlap"] = map["dst_global_y_size_overlap"]

            self.map_mpi[slice_index] = map

            slice_index = slice_index + 1

    def update_mpi_map(self):
        """
        Update the MPI map for parallel processing.
            This method recalculates the source and overlap slices for each process based on the current grid.

        Examples
        --------
        >>> coverage.update_mpi_map()
        """
        if self.is_regular_grid(type="source"):

            idx = np.where((self.source_global_axis_x >= np.min(self.read_axis_x(type="target", with_overlap=False))) &
                           (self.source_global_axis_x <= np.max(self.read_axis_x(type="target", with_overlap=False))))

            xmin = np.min(idx[0])
            xmax = np.max(idx[0]) + 1

            idx = np.where((self.source_global_axis_y >= np.min(self.read_axis_y(type="target", with_overlap=False))) &
                           (self.source_global_axis_y <= np.max(self.read_axis_y(type="target", with_overlap=False))))

            ymin = np.min(idx[0])
            ymax = np.max(idx[0]) + 1

        else:

            idx = np.where(
                (self.source_global_axis_x >= np.min(self.read_axis_x(type="target", with_overlap=False))) &
                (self.source_global_axis_x <= np.max(self.read_axis_x(type="target", with_overlap=False))) &
                (self.source_global_axis_y >= np.min(self.read_axis_y(type="target", with_overlap=False))) &
                (self.source_global_axis_y <= np.max(self.read_axis_y(type="target", with_overlap=False))))

            ymin = np.min(idx[0])
            ymax = np.max(idx[0]) + 1
            xmin = np.min(idx[1])
            xmax = np.max(idx[1]) + 1

        # Version 2
        # SRC GLOBAL
        self.map_mpi[self.rank]["src_global_x"] = np.s_[xmin:xmax]
        self.map_mpi[self.rank]["src_global_x_size"] = xmax - xmin
        self.map_mpi[self.rank]["src_global_y"] = np.s_[ymin:ymax]
        self.map_mpi[self.rank]["src_global_y_size"] = ymax - ymin

        dst_global_x_min_overlap = max(0, self.map_mpi[self.rank][
            "src_global_x"].start - Coverage.HORIZONTAL_OVERLAPING_SIZE)
        dst_global_x_max_overlap = min(self.source_global_x_size,
                                       self.map_mpi[self.rank][
                                           "src_global_x"].stop + Coverage.HORIZONTAL_OVERLAPING_SIZE)
        self.map_mpi[self.rank]["src_global_x_overlap"] = np.s_[
                                                          dst_global_x_min_overlap:dst_global_x_max_overlap]

        dst_global_y_min_overlap = max(0, self.map_mpi[self.rank][
            "src_global_y"].start - Coverage.HORIZONTAL_OVERLAPING_SIZE)
        dst_global_y_max_overlap = min(self.source_global_y_size,
                                       self.map_mpi[self.rank][
                                           "src_global_y"].stop + Coverage.HORIZONTAL_OVERLAPING_SIZE)
        self.map_mpi[self.rank]["src_global_y_overlap"] = np.s_[
                                                          dst_global_y_min_overlap:dst_global_y_max_overlap]

        self.map_mpi[self.rank]["src_global_x_size_overlap"] = self.map_mpi[self.rank][
                                                                   "src_global_x_overlap"].stop - \
                                                               self.map_mpi[self.rank][
                                                                   "src_global_x_overlap"].start
        self.map_mpi[self.rank]["src_global_y_size_overlap"] = self.map_mpi[self.rank][
                                                                   "src_global_y_overlap"].stop - \
                                                               self.map_mpi[self.rank][
                                                                   "src_global_y_overlap"].start

        self.map_mpi[self.rank]["src_local_x_size"] = xmax - xmin
        self.map_mpi[self.rank]["src_local_y_size"] = ymax - ymin
        self.map_mpi[self.rank]["src_local_x"] = np.s_[0:self.map_mpi[self.rank]["src_local_x_size"]]
        self.map_mpi[self.rank]["src_local_y"] = np.s_[0:self.map_mpi[self.rank]["src_local_y_size"]]

        # OVERLAP
        self.map_mpi[self.rank]["src_local_x_size_overlap"] = self.map_mpi[self.rank][
            "src_global_x_size_overlap"]
        self.map_mpi[self.rank]["src_local_y_size_overlap"] = self.map_mpi[self.rank][
            "src_global_y_size_overlap"]

        self.map_mpi[self.rank]["src_local_x_overlap"] = np.s_[
                                                         0:self.map_mpi[self.rank]["src_local_x_size_overlap"]]
        self.map_mpi[self.rank]["src_local_y_overlap"] = np.s_[
                                                         0:self.map_mpi[self.rank]["src_local_y_size_overlap"]]

    def create_threading_map(self):
        """
        Create the multi threading map for parallel processing.

        Examples
        --------
        >>> coverage.create_threading_map()
        """
        self.threading_map = np.empty([self.size, self.threads_number], dtype=object)
        target_sample = (self.map_mpi[self.rank]["dst_local_y_size"], self.map_mpi[self.rank]["dst_local_x_size"])

        # Split the axes
        target_slices = shape_split(target_sample, self.threads_number, axis=[0, 0])

        # If we can divide the grid dimensions with the number of threads,
        # we set the threads number by the number of slice
        self.threads_number =len(target_slices.flatten())

        slice_index = 0
        for slyce in target_slices.flatten():
            slice = tuple(slyce)

            map = {}

            #### Destination grid ###
            map["dst_global_x"] = slice[1]
            map["dst_global_y"] = slice[0]

            map["dst_local_x_size"] = map["dst_global_x"].stop - map["dst_global_x"].start
            map["dst_local_y_size"] = map["dst_global_y"].stop - map["dst_global_y"].start

            # Compute overlap
            dst_global_x_min_overlap = max(0, map["dst_global_x"].start - Coverage.HORIZONTAL_OVERLAPING_SIZE)
            dst_global_x_max_overlap = min(self.map_mpi[self.rank]["dst_local_x_size"],
                                           map["dst_global_x"].stop + Coverage.HORIZONTAL_OVERLAPING_SIZE)
            map["dst_global_x_overlap"] = np.s_[dst_global_x_min_overlap:dst_global_x_max_overlap]

            dst_global_y_min_overlap = max(0, map["dst_global_y"].start - Coverage.HORIZONTAL_OVERLAPING_SIZE)
            dst_global_y_max_overlap = min(self.map_mpi[self.rank]["dst_local_y_size"],
                                           map["dst_global_y"].stop + Coverage.HORIZONTAL_OVERLAPING_SIZE)
            map["dst_global_y_overlap"] = np.s_[dst_global_y_min_overlap:dst_global_y_max_overlap]

            map["dst_global_x_size_overlap"] = map["dst_global_x_overlap"].stop - map["dst_global_x_overlap"].start
            map["dst_global_y_size_overlap"] = map["dst_global_y_overlap"].stop - map["dst_global_y_overlap"].start

            # Compute dest local grid
            dst_local_x_min = map["dst_global_x"].start
            dst_local_x_max = map["dst_global_x"].stop
            dst_local_y_min = map["dst_global_y"].start
            dst_local_y_max = map["dst_global_y"].stop

            map["dst_local_x"] = np.s_[dst_local_x_min:dst_local_x_max]
            map["dst_local_y"] = np.s_[dst_local_y_min:dst_local_y_max]

            ### Source grille ###
            map["src_global_x"] = map["dst_global_x"]
            map["src_global_y"] = map["dst_global_y"]

            map["src_global_x_overlap"] = map["dst_global_x_overlap"]
            map["src_global_y_overlap"] = map["dst_global_y_overlap"]

            map["src_local_x"] = map["dst_local_x"]
            map["src_local_y"] = map["dst_local_y"]

            map["src_local_x_size"] = map["dst_local_x_size"]
            map["src_local_y_size"] = map["dst_local_y_size"]

            map["src_local_x_size_overlap"] = map["dst_global_x_size_overlap"]
            map["src_local_y_size_overlap"] = map["dst_global_y_size_overlap"]

            self.threading_map[self.rank][slice_index] = map

            slice_index = slice_index + 1

    def update_threading_map(self):
        """
        Update source grid for threading map

        Examples
        --------
        >>> coverage.update_mpi_map()
        """
        if self.is_regular_grid(type="source"):

            idx = np.where((self.source_global_axis_x >= np.min(self.read_axis_x(type="target", with_overlap=False))) &
                           (self.source_global_axis_x <= np.max(self.read_axis_x(type="target", with_overlap=False))))

            xmin = np.min(idx[0])
            xmax = np.max(idx[0]) + 1

            idx = np.where((self.source_global_axis_y >= np.min(self.read_axis_y(type="target", with_overlap=False))) &
                           (self.source_global_axis_y <= np.max(self.read_axis_y(type="target", with_overlap=False))))

            ymin = np.min(idx[0])
            ymax = np.max(idx[0]) + 1

        else:

            idx = np.where(
                (self.source_global_axis_x >= np.min(self.read_axis_x(type="target", with_overlap=False))) &
                (self.source_global_axis_x <= np.max(self.read_axis_x(type="target", with_overlap=False))) &
                (self.source_global_axis_y >= np.min(self.read_axis_y(type="target", with_overlap=False))) &
                (self.source_global_axis_y <= np.max(self.read_axis_y(type="target", with_overlap=False))))

            ymin = np.min(idx[0])
            ymax = np.max(idx[0]) + 1
            xmin = np.min(idx[1])
            xmax = np.max(idx[1]) + 1

        for thread in range(0,self.threads_number):
            # Version 2
            # SRC GLOBAL
            self.threading_map[self.rank][thread]["src_global_x"] = np.s_[xmin:xmax]
            self.threading_map[self.rank][thread]["src_global_x_size"] = xmax - xmin
            self.threading_map[self.rank][thread]["src_global_y"] = np.s_[ymin:ymax]
            self.threading_map[self.rank][thread]["src_global_y_size"] = ymax - ymin

            dst_global_x_min_overlap = max(0, self.threading_map[self.rank][thread][
                "src_global_x"].start - Coverage.HORIZONTAL_OVERLAPING_SIZE)
            dst_global_x_max_overlap = min(self.map_mpi[self.rank]["dst_local_x_size"],
                                           self.threading_map[self.rank][thread][
                                               "src_global_x"].stop + Coverage.HORIZONTAL_OVERLAPING_SIZE)
            self.threading_map[self.rank][thread]["src_global_x_overlap"] = np.s_[
                                                              dst_global_x_min_overlap:dst_global_x_max_overlap]

            dst_global_y_min_overlap = max(0, self.threading_map[self.rank][thread][
                "src_global_y"].start - Coverage.HORIZONTAL_OVERLAPING_SIZE)
            dst_global_y_max_overlap = min(self.map_mpi[self.rank]["dst_local_y_size"],
                                           self.threading_map[self.rank][thread][
                                               "src_global_y"].stop + Coverage.HORIZONTAL_OVERLAPING_SIZE)
            self.threading_map[self.rank][thread]["src_global_y_overlap"] = np.s_[
                                                              dst_global_y_min_overlap:dst_global_y_max_overlap]

            self.threading_map[self.rank][thread]["src_global_x_size_overlap"] = self.threading_map[self.rank][thread][
                                                                       "src_global_x_overlap"].stop - \
                                                                   self.threading_map[self.rank][thread][
                                                                       "src_global_x_overlap"].start
            self.threading_map[self.rank][thread]["src_global_y_size_overlap"] = self.threading_map[self.rank][thread][
                                                                       "src_global_y_overlap"].stop - \
                                                                   self.threading_map[self.rank][thread][
                                                                       "src_global_y_overlap"].start

            self.threading_map[self.rank][thread]["src_local_x_size"] = xmax - xmin
            self.threading_map[self.rank][thread]["src_local_y_size"] = ymax - ymin
            self.threading_map[self.rank][thread]["src_local_x"] = np.s_[0:self.threading_map[self.rank][thread]["src_local_x_size"]]
            self.threading_map[self.rank][thread]["src_local_y"] = np.s_[0:self.threading_map[self.rank][thread]["src_local_y_size"]]

            # OVERLAP
            self.threading_map[self.rank][thread]["src_local_x_size_overlap"] = self.threading_map[self.rank][thread][
                "src_global_x_size_overlap"]
            self.threading_map[self.rank][thread]["src_local_y_size_overlap"] = self.threading_map[self.rank][thread][
                "src_global_y_size_overlap"]

            self.threading_map[self.rank][thread]["src_local_x_overlap"] = np.s_[
                                                             0:self.threading_map[self.rank][thread]["src_local_x_size_overlap"]]
            self.threading_map[self.rank][thread]["src_local_y_overlap"] = np.s_[
                                                             0:self.threading_map[self.rank][thread]["src_local_y_size_overlap"]]

    # Read metadata
    def read_metadata(self):
        """
        Read the metadata from the file if the reader contains a read_metadata() function.
        Examples
        --------
        >>> coverage.read_metadata()
        """
        if "read_metadata" in dir(self.reader):
            m = self.reader.read_metadata()

    def get_x_size(self, type="target", with_overlap=False):
        """
        Get the x size.

        Parameters
        ----------
        type : str, optional
            Type of the grid ("target", "source", "target_global", "source_global").
        with_overlap : bool, optional
            Whether to include overlap.

        Returns
        -------
        int
            The x size.
        Examples
        --------
        >>> coverage.get_x_size()
        120
        >>> coverage.get_x_size(type="source", with_overlap=True)
        124
        """
        if type == "target_global":
            return self.target_global_x_size
        elif type == "source_global":
            return self.source_global_x_size
        elif type == "source" and with_overlap is True:
            return self.map_mpi[self.rank]["src_local_x_size_overlap"]
        elif type == "source" and with_overlap is False:
            return self.map_mpi[self.rank]["src_local_x_size"]
        elif type == "target" and with_overlap is True:
            return self.map_mpi[self.rank]["dst_local_x_size_overlap"]
        else:
            return self.map_mpi[self.rank]["dst_local_x_size"]

    def get_y_size(self, type="target", with_overlap=False):
        """
        Get the y size.

        Parameters
        ----------
        type : str, optional
            Type of the grid ("target", "source", "target_global", "source_global").
        with_overlap : bool, optional
            Whether to include overlap.

        Returns
        -------
        int
            The y size.

        Examples
        --------
        >>> coverage.get_y_size()
        80
        >>> coverage.get_y_size(type="source", with_overlap=True)
        84
        """
        if type == "target_global":
            return self.target_global_y_size
        elif type == "source_global":
            return self.source_global_y_size
        elif type == "source" and with_overlap is True:
            return self.map_mpi[self.rank]["src_local_y_size_overlap"]
        elif type == "source" and with_overlap is False:
            return self.map_mpi[self.rank]["src_local_y_size"]
        elif type == "target" and with_overlap is True:
            return self.map_mpi[self.rank]["dst_local_y_size_overlap"]
        else:
            return self.map_mpi[self.rank]["dst_local_y_size"]

    def is_regular_grid(self, type="target"):
        """
        Check if the grid is regular.

        Parameters
        ----------
        type : str, optional
            Type of the grid ("target", "source").

        Returns
        -------
        bool
            True if the grid is regular, False otherwise.
        Examples
        --------
        >>> coverage.is_regular_grid()
        True
        >>> coverage.is_regular_grid(type="source")
        False
        """
        if type == "target":
            return self.target_regular_grid
        else:
            return self.source_regular_grid

    # Axis        
    def read_axis_x(self, type="target", with_overlap=False):
        """
        Return the values (often longitude) of the x axis.

        Parameters
        ----------
        type : str, optional
            Type of the grid ("target", "source", "target_global", "source_global").
        with_overlap : bool, optional
            Whether to include overlap.

        Returns
        -------
        array
            A one- or two-dimensional array of x axis values (often longitude): [x] or [y, x], depending on the grid type.

        Examples
        --------
        >>> x_axis = coverage.read_axis_x()
        >>> x_axis = coverage.read_axis_x(type="source", with_overlap=True)
        """

        if type == "target_global":
            return self.target_global_axis_x

        elif type == "source_global":
            return self.source_global_axis_x

        elif type == "source" and with_overlap is True:
            return self.reader.read_axis_x(self.map_mpi[self.rank]["src_global_x_overlap"].start,
                                           self.map_mpi[self.rank]["src_global_x_overlap"].stop,
                                           self.map_mpi[self.rank]["src_global_y_overlap"].start,
                                           self.map_mpi[self.rank]["src_global_y_overlap"].stop)
        elif type == "source" and with_overlap is False:
            return self.reader.read_axis_x(self.map_mpi[self.rank]["src_global_x"].start,
                                           self.map_mpi[self.rank]["src_global_x"].stop,
                                           self.map_mpi[self.rank]["src_global_y"].start,
                                           self.map_mpi[self.rank]["src_global_y"].stop)

        elif type == "target" and with_overlap is True:

            if self.is_regular_grid():
                return self.target_global_axis_x[self.map_mpi[self.rank]["dst_global_x_overlap"]]
            else:
                return self.target_global_axis_x[self.map_mpi[self.rank]["dst_global_y_overlap"],
                self.map_mpi[self.rank]["dst_global_x_overlap"]]
        else:

            if self.is_regular_grid():
                return self.target_global_axis_x[self.map_mpi[self.rank]["dst_global_x"]]
            else:
                return self.target_global_axis_x[
                    self.map_mpi[self.rank]["dst_global_y"], self.map_mpi[self.rank]["dst_global_x"]]

    def read_axis_y(self, type="target", with_overlap=False):
        """
        Return the values (often latitude) of the y axis.

        Parameters
        ----------
        type : str, optional
            Type of the grid ("target", "source", "target_global", "source_global").
        with_overlap : bool, optional
            Whether to include overlap.

        Returns
        -------
        array
            A one- or two-dimensional array of y axis values (often latitude): [x] or [y, x], depending on the grid type.

        Examples
        --------
        >>> y_axis = coverage.read_axis_y()
        >>> y_axis = coverage.read_axis_y(type="source", with_overlap=True)
        """
        if type == "target_global":
            return self.target_global_axis_y

        elif type == "source_global":
            return self.source_global_axis_y

        elif type == "source" and with_overlap is True:
            return self.reader.read_axis_y(self.map_mpi[self.rank]["src_global_x_overlap"].start,
                                           self.map_mpi[self.rank]["src_global_x_overlap"].stop,
                                           self.map_mpi[self.rank]["src_global_y_overlap"].start,
                                           self.map_mpi[self.rank]["src_global_y_overlap"].stop)
        elif type == "source" and with_overlap is False:
            return self.reader.read_axis_y(self.map_mpi[self.rank]["src_global_x"].start,
                                           self.map_mpi[self.rank]["src_global_x"].stop,
                                           self.map_mpi[self.rank]["src_global_y"].start,
                                           self.map_mpi[self.rank]["src_global_y"].stop)

        elif type == "target" and with_overlap is True:

            if self.is_regular_grid():
                return self.target_global_axis_y[self.map_mpi[self.rank]["dst_global_y_overlap"]]
            else:
                return self.target_global_axis_y[self.map_mpi[self.rank]["dst_global_y_overlap"],
                self.map_mpi[self.rank]["dst_global_x_overlap"]]
        else:

            if self.is_regular_grid():
                return self.target_global_axis_y[self.map_mpi[self.rank]["dst_global_y"]]
            else:
                return self.target_global_axis_y[
                    self.map_mpi[self.rank]["dst_global_y"], self.map_mpi[self.rank]["dst_global_x"]]

    def find_point_index(self, target_lon, target_lat, decimal_tolerance=5, method="classic", only_mask_value=True,
                         type="source"):
        """
                Find the index of the closest point to the given point.

                Parameters
                ----------
                target_lon : float
                    Longitude of the target point.
                target_lat : float
                    Latitude of the target point.
                decimal_tolerance : int, optional
                    Tolerance for rounding coordinates.
                method : str, optional
                    Method for calculation ("classic", "quick").
                only_mask_value : bool, optional
                    Whether to consider only mask values.
                type : str, optional
                    Type of the grid ("source", "source_global").

                Returns
                -------
                list
                    A list containing:
                    [0] : the x index of the closest point
                    [1] : the y index of the closest point
                    [2] : the longitude of the closest point
                    [3] : the latitude of the closest point
                    [4] : the distance of the closest point in kilometers.

                Examples
                --------
                >>> idx = coverage.find_point_index(5.0, 45.0)
                >>> idx = coverage.find_point_index(5.0, 45.0, method="quick", type="source_global")
                """
        lon = self.read_axis_x(type="source", with_overlap=False)
        lat = self.read_axis_y(type="source", with_overlap=False)

        if self.check_point_is_inside(target_lon, target_lat, lon, lat, tolerance=decimal_tolerance):
            try:
                mask = self.read_variable_2D_sea_binary_mask(type="source", with_overlap=False)
            except NotImplementedError:
                logging.warning("No 2D sea binary mask found")
                # mask = np.ones([self.source_global_y_size, self.source_global_x_size])
                only_mask_value = False

            dist = np.zeros([self.get_y_size(type="source", with_overlap=False),
                             self.get_x_size(type="source", with_overlap=False)])
            dist[:] = 10000000

            if method == "classic":
                for x in range(0, self.get_x_size(type="source", with_overlap=False)):
                    for y in range(0, self.get_y_size(type="source", with_overlap=False)):

                        if only_mask_value:
                            if (mask[y, x] == 1):  # =Terre
                                if self.reader.is_regular_grid():
                                    dist[y, x] = distance_on_unit_sphere(target_lon, target_lat, lon[x], lat[y])
                                else:
                                    dist[y, x] = distance_on_unit_sphere(target_lon, target_lat, lon[y, x], lat[y, x])
                        else:
                            if self.reader.is_regular_grid():
                                dist[y, x] = distance_on_unit_sphere(target_lon, target_lat, lon[x], lat[y])
                            else:
                                dist[y, x] = distance_on_unit_sphere(target_lon, target_lat, lon[y, x], lat[y, x])

                nearest_y_index, nearest_x_index = np.where(dist == np.min(dist))

                if len(nearest_y_index) == 0 or len(nearest_x_index) == 0:
                    logging.error("no neaarest point found")
                    raise RuntimeError("No nearest point found")

                nearest_x_index = nearest_x_index[0]
                nearest_y_index = nearest_y_index[0]
                min_dist = dist[nearest_y_index, nearest_x_index]

                if self.is_regular_grid(type="source"):
                    nearest_lon = lon[nearest_x_index]
                    nearest_lat = lat[nearest_y_index]
                else:
                    nearest_lon = lon[nearest_y_index, nearest_x_index]
                    nearest_lat = lat[nearest_y_index, nearest_x_index]

                if type == "source":
                    return [nearest_x_index, nearest_y_index, nearest_lon, nearest_lat, min_dist]
                elif type == "source_global":
                    return [self.map_mpi[self.rank]["src_global_x"].start + nearest_x_index,
                            self.map_mpi[self.rank]["src_global_y"].start + nearest_y_index, nearest_lon, nearest_lat,
                            min_dist]
                else:
                    raise ValueError("Type doesn't match [source, source_global]")

            elif method == "quick":

                if self.is_regular_grid():

                    # Longitude : on cherche l'index le plus proche
                    array = np.asarray(lon)
                    nearest_x_index = (np.abs(array - target_lon)).argmin()

                    # Latitude : on cherche l'index le plus proche
                    array = np.asarray(lat)
                    nearest_y_index = (np.abs(array - target_lat)).argmin()

                    nearest_lon = lon[nearest_x_index]
                    nearest_lat = lat[nearest_y_index]

                    min_dist = distance_on_unit_sphere(target_lon, target_lat, nearest_lon, nearest_lat)

                    if type == "source":
                        return [nearest_x_index, nearest_y_index, nearest_lon, nearest_lat, min_dist]
                    elif type == "source_global":
                        return [self.map_mpi[self.rank]["src_global_x"].start + nearest_x_index,
                                self.map_mpi[self.rank]["src_global_y"].start + nearest_y_index, nearest_lon,
                                nearest_lat, min_dist]
                    else:
                        raise ValueError("Type doesn't match [source, source_global]")

                else:
                    raise NotImplementedError("Method " + str(method) + " is not implemented for regular grid.")

            else:
                raise RuntimeError("Method " + str(method) + " is not implemented yet.")
        else:
            logging.warning("Point is outside the rank n°" + str(self.rank))
            raise NotFoundInRankError(self.rank, "Point is outside the rank")

    # Variables
    #################
    # HYDRO
    # 2D
    #################
    def read_variable_bathymetry(self,current_thread:int=0):
        """
        Read the bathymetry over the entire coverage.

        Returns
        -------
        array
            A two-dimensional array [y, x].

        Examples
        --------
        >>> bathymetry = coverage.read_variable_bathymetry()
        """
        data = self.reader.read_variable_bathymetry(
            self.threading_map[self.rank][current_thread]["src_global_x_overlap"].start,
            self.threading_map[self.rank][current_thread]["src_global_x_overlap"].stop,
            self.threading_map[self.rank][current_thread]["src_global_y_overlap"].start,
            self.threading_map[self.rank][current_thread]["src_global_y_overlap"].stop)

        if self.horizontal_resampling:
            data = resample_2d_to_grid(self.read_axis_x(type="source", with_overlap=True),
                                       self.read_axis_y(type="source", with_overlap=True),
                                       self.read_axis_x(type="target", with_overlap=True),
                                       self.read_axis_y(type="target", with_overlap=True),
                                       data,
                                       Coverage.HORIZONTAL_INTERPOLATION_METHOD)

        return data[self.threading_map[self.rank][current_thread]["dst_local_y"], self.threading_map[self.rank][current_thread]["dst_local_x"]]

    def read_variable_topography(self):
        """
        Read the topography over the entire coverage.

        Returns
        -------
        array
            A two-dimensional array [y, x].

        Examples
        --------
        >>> topography = coverage.read_variable_topography()
        """
        data = self.reader.read_variable_topography(
            self.map_mpi[self.rank]["src_global_x_overlap"].start,
            self.map_mpi[self.rank]["src_global_x_overlap"].stop,
            self.map_mpi[self.rank]["src_global_y_overlap"].start,
            self.map_mpi[self.rank]["src_global_y_overlap"].stop)

        if self.horizontal_resampling:
            data = resample_2d_to_grid(self.read_axis_x(type="source", with_overlap=True),
                                       self.read_axis_y(type="source", with_overlap=True),
                                       self.read_axis_x(type="target", with_overlap=True),
                                       self.read_axis_y(type="target", with_overlap=True),
                                       data,
                                       Coverage.HORIZONTAL_INTERPOLATION_METHOD)

        return data[self.map_mpi[self.rank]["dst_local_y"], self.map_mpi[self.rank]["dst_local_x"]]

    def read_variable_mesh_size(self):
        """
        Read the mesh size over the entire coverage.

        Returns
        -------
        array
            A two-dimensional array [y, x].

        Examples
        --------
        >>> mesh_size = coverage.read_variable_mesh_size()
        """
        data = self.reader.read_variable_mesh_size(
            self.map_mpi[self.rank]["src_global_x_overlap"].start,
            self.map_mpi[self.rank]["src_global_x_overlap"].stop,
            self.map_mpi[self.rank]["src_global_y_overlap"].start,
            self.map_mpi[self.rank]["src_global_y_overlap"].stop)

        if self.horizontal_resampling:
            data = resample_2d_to_grid(self.read_axis_x(type="source", with_overlap=True),
                                       self.read_axis_y(type="source", with_overlap=True),
                                       self.read_axis_x(type="target", with_overlap=True),
                                       self.read_axis_y(type="target", with_overlap=True),
                                       data,
                                       Coverage.HORIZONTAL_INTERPOLATION_METHOD)

        return data[self.map_mpi[self.rank]["dst_local_y"], self.map_mpi[self.rank]["dst_local_x"]]

    def read_variable_x_mesh_size(self):
        """
        Read the mesh size on the x axis.

        Returns
        -------
        array
            A two-dimensional array [y, x].

        Examples
        --------
        >>> x_mesh_size = coverage.read_variable_x_mesh_size()
        """
        data = self.reader.read_variable_x_mesh_size(
            self.map_mpi[self.rank]["src_global_x_overlap"].start,
            self.map_mpi[self.rank]["src_global_x_overlap"].stop,
            self.map_mpi[self.rank]["src_global_y_overlap"].start,
            self.map_mpi[self.rank]["src_global_y_overlap"].stop)

        if self.horizontal_resampling:
            data = resample_2d_to_grid(self.read_axis_x(type="source", with_overlap=True),
                                       self.read_axis_y(type="source", with_overlap=True),
                                       self.read_axis_x(type="target", with_overlap=True),
                                       self.read_axis_y(type="target", with_overlap=True),
                                       data,
                                       Coverage.HORIZONTAL_INTERPOLATION_METHOD)

        return data[self.map_mpi[self.rank]["dst_local_y"], self.map_mpi[self.rank]["dst_local_x"]]

    def read_variable_y_mesh_size(self):
        """
        Read the mesh size on the y axis.

        Returns
        -------
        array
            A two-dimensional array [y, x].

        Examples
        --------
        >>> y_mesh_size = coverage.read_variable_y_mesh_size()
        """
        data = self.reader.read_variable_y_mesh_size(
            self.map_mpi[self.rank]["src_global_x_overlap"].start,
            self.map_mpi[self.rank]["src_global_x_overlap"].stop,
            self.map_mpi[self.rank]["src_global_y_overlap"].start,
            self.map_mpi[self.rank]["src_global_y_overlap"].stop)

        if self.horizontal_resampling:
            data = resample_2d_to_grid(self.read_axis_x(type="source", with_overlap=True),
                                       self.read_axis_y(type="source", with_overlap=True),
                                       self.read_axis_x(type="target", with_overlap=True),
                                       self.read_axis_y(type="target", with_overlap=True),
                                       data,
                                       Coverage.HORIZONTAL_INTERPOLATION_METHOD)

        return data[self.map_mpi[self.rank]["dst_local_y"], self.map_mpi[self.rank]["dst_local_x"]]

    def read_variable_2D_sea_binary_mask(self, type="target", with_overlap=False):
        """
        Read the land/sea mask over the entire coverage.

        Parameters
        ----------
        type : str, optional
            Type of the grid ("target", "source").
        with_overlap : bool, optional
            Whether to include overlap.

        Returns
        -------
        array
            A two-dimensional array [y, x].
            0 = Land
            1 = Sea
        Examples
        --------
        >>> mask = coverage.read_variable_2D_sea_binary_mask()
        """
        data = self.reader.read_variable_2D_sea_binary_mask(
            self.map_mpi[self.rank]["src_global_x_overlap"].start,
            self.map_mpi[self.rank]["src_global_x_overlap"].stop,
            self.map_mpi[self.rank]["src_global_y_overlap"].start,
            self.map_mpi[self.rank]["src_global_y_overlap"].stop)

        if type == "source" and with_overlap is True:
            return data[self.map_mpi[self.rank]["src_local_y_overlap"], self.map_mpi[self.rank]["src_local_x_overlap"]]
        elif type == "source" and with_overlap is False:
            return data[self.map_mpi[self.rank]["src_local_y"], self.map_mpi[self.rank]["src_local_x"]]
        elif type == "target":

            if self.horizontal_resampling:
                data = resample_2d_to_grid(self.read_axis_x(type="source", with_overlap=True),
                                           self.read_axis_y(type="source", with_overlap=True),
                                           self.read_axis_x(type="target", with_overlap=True),
                                           self.read_axis_y(type="target", with_overlap=True),
                                           data,
                                           Coverage.HORIZONTAL_INTERPOLATION_METHOD)

            return data[self.map_mpi[self.rank]["dst_local_y"], self.map_mpi[self.rank]["dst_local_x"]]

    def read_variable_Ha(self):
        """
        Read the amplitude of the reanalysis.

        Returns
        -------
        array
            A two-dimensional array [y, x].

        Examples
        --------
        >>> Ha = coverage.read_variable_Ha()
        """
        data = self.reader.read_variable_Ha(
            self.map_mpi[self.rank]["src_global_x_overlap"].start,
            self.map_mpi[self.rank]["src_global_x_overlap"].stop,
            self.map_mpi[self.rank]["src_global_y_overlap"].start,
            self.map_mpi[self.rank]["src_global_y_overlap"].stop)

        if self.horizontal_resampling:
            data = resample_2d_to_grid(self.read_axis_x(type="source", with_overlap=True),
                                       self.read_axis_y(type="source", with_overlap=True),
                                       self.read_axis_x(type="target", with_overlap=True),
                                       self.read_axis_y(type="target", with_overlap=True),
                                       data,
                                       Coverage.HORIZONTAL_INTERPOLATION_METHOD)

        return data[self.map_mpi[self.rank]["dst_local_y"], self.map_mpi[self.rank]["dst_local_x"]]
