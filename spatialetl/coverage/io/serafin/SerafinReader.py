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
"""!
Read/Write Serafin files and manipulate associated data.

Handles Serafin file with:
  - single and double precision
  - big and little endian

The sizes (file, header, frame) are in bytes (8 bits).
"""
from __future__ import division, print_function, absolute_import
from spatialetl.coverage.io.CoverageReader import CoverageReader
import copy
import numpy as np
import os
import struct
import logging
from datetime import datetime,timedelta
from spatialetl.coverage.io.serafin.SerafinHeader import SerafinHeader

class SerafinReader(CoverageReader):

    def __init__(self, myFilename,language="fr"):
        CoverageReader.__init__(self, myFilename);
        self.language = language

        self.file = open(self.filename, "rb")
        self.file_size =os.path.getsize(self.filename)
        self.header = SerafinHeader(lang=self.language)
        self.header.from_file(self.file, self.file_size)

        self.x_size = np.shape(self.header.x)[0]
        self.y_size = np.shape(self.header.y)[0]

        self.time_ref=datetime(self.header.date[0],self.header.date[1],self.header.date[2],self.header.date[3],self.header.date[4],self.header.date[5])

    def close(self):
        self.file.close()

    def read_axis_x(self):
        return self.header.x

    def read_axis_y(self):
        return self.header.y

    def read_axis_t(self, timestamp):
        result = []

        if self.header is None:
            raise ValueError('Cannot read time without any header (forgot read_header ?)')
        logging.debug('Reading the time series from the file')
        self.file.seek(self.header.header_size, 0)

        for _ in range(self.header.nb_frames):
            self.file.read(4)
            d = timedelta(seconds=self.header.unpack_float(self.file.read(self.header.float_size), 1)[0])
            result.append(self.time_ref+d)
            self.file.read(4)
            self.file.seek(self.header.frame_size - 8 - self.header.float_size, 1)

        if timestamp == 1:
            return [(t - TimeCoverage.TIME_DATUM).total_seconds() \
                    for t in result];
        else:
            return result;


    # Variables
    def read_variable_longitude(self):
        return self.read_axis_x()

    def read_variable_latitude(self):
        return self.read_axis_y()

    def read_variable_time(self):
        return self.read_axis_t(timestamp=0)

    def _get_var_index(self, var_ID):
        """!
        @brief Handle data request by variable ID
        @param var_ID <str>: the ID of the requested variable
        @return index <int> the index of the frame (0-based)
        """
        if self.header is None:
            raise ValueError('Cannot extract variable from empty list (forgot read_header ?)')
        try:
            index = self.header.var_IDs.index(var_ID)
        except ValueError:
            raise ValueError('Variable ID %s not found' % var_ID)
        return index

    def read_variable_bathymetry(self):
        time_index=0
        var_ID="WATER DEPTH"
        pos_var = self._get_var_index(var_ID)
        self.file.seek(self.header.header_size + time_index * self.header.frame_size
                       + 8 + self.header.float_size + pos_var * (8 + self.header.float_size * self.header.nb_nodes), 0)
        self.file.read(4)
        data = np.array(self.header.unpack_float(self.file.read(self.header.float_size * self.header.nb_nodes),
                                                 self.header.nb_nodes), dtype=self.header.np_float_type)

        return data


    def read_var_in_frame(self, time_index, var_ID):
        """!
        @brief Read a single variable in a frame
        @param time_index <int>: the index of the frame (0-based)
        @param var_ID <str>: variable ID
        @return <numpy 1D-array>: values of the variables, of length equal to the number of nodes
        """
        if time_index < 0:
            raise SerafinRequestError('Impossible to read a negative time index!')
        logger.debug('Reading variable %s at frame %i' % (var_ID, time_index))
        pos_var = self._get_var_index(var_ID)
        self.file.seek(self.header.header_size + time_index * self.header.frame_size
                       + 8 + self.header.float_size + pos_var * (8 + self.header.float_size * self.header.nb_nodes), 0)
        self.file.read(4)
        return np.array(self.header.unpack_float(self.file.read(self.header.float_size * self.header.nb_nodes),
                                                 self.header.nb_nodes), dtype=self.header.np_float_type)

    def read_var_in_frame_as_3d(self, time_index, var_ID):
        """!
        @brief Read a single variable in a 3D frame
        @param time_index <int>: the index of the frame (0-based)
        @param var_ID <str>: variable ID
        @return <numpy 2D-array>: values of the variables with shape (planes number, number of 2D nodes)
        """
        if self.header.is_2d:
            raise ValueError('Reading values as 3D is only possible in 3D!')
        new_shape = (self.header.nb_planes, self.header.nb_nodes_2d)
        return self.read_var_in_frame(time_index, var_ID).reshape(new_shape)

    def read_var_in_frame_at_layer(self, time_index, var_ID, iplan):
        """!
        @brief Read a single variable in a frame at specific layer
        @param time_index <int>: the index of the frame (0-based)
        @param var_ID <str>: variable ID
        @param iplan <int>: 1-based index of layer
        @return <numpy 1D-array>: values of the variables, of length equal to the number of nodes
        """
        if self.header.is_2d:
            raise ValueError('Extracting values at a specific layer is only possible in 3D!')
        if iplan < 1 or iplan > self.header.nb_planes:
            raise SerafinRequestError('Layer %i is not inside [1, %i]' % (iplan, self.header.nb_planes))
        return self.read_var_in_frame_as_3d(time_index, var_ID)[iplan + 1]