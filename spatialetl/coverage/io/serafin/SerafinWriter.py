# #! /usr/bin/env python2.7
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
from spatialetl.coverage.io.CoverageWriter import CoverageWriter

class SerafinWriter(CoverageWriter):

    def __init__(self, cov,myFile,depths,overwrite=False):
        CoverageWriter.__init__(self,cov,myFile);

        mode = 'wb' if overwrite else 'xb'
        self.language = language

        logger.info('Writing the output file: "%s"' % filename)

        logger.debug('Writing header with %i nodes, %i elements %s and %i variable %s' %
                     (header.nb_nodes, header.nb_elements,
                      '' if header.is_2d else ', %i layers' % header.nb_planes,
                      header.nb_var, ['', 's'][header.nb_var > 1]))

        # Title and file type
        self.file.write(header.pack_int(80))
        self.file.write(header.title)
        self.file.write(header.file_format)
        self.file.write(header.pack_int(80))

        # Number of variables
        self.file.write(header.pack_int(8))
        self.file.write(header.pack_int(header.nb_var))
        self.file.write(header.pack_int(header.nb_var_quadratic))
        self.file.write(header.pack_int(8))

        # Variable names and units
        for j in range(header.nb_var):
            self.file.write(header.pack_int(2 * 16))
            self.file.write(header.var_names[j].ljust(16))
            self.file.write(header.var_units[j].ljust(16))
            self.file.write(header.pack_int(2 * 16))

        # Date
        self.file.write(header.pack_int(10 * 4))
        self.file.write(header.pack_int(*header.params, nb=10))
        self.file.write(header.pack_int(10 * 4))
        if header.params[-1] == 1:
            self.file.write(header.pack_int(6 * 4))
            self.file.write(header.pack_int(*header.date, nb=6))
            self.file.write(header.pack_int(6 * 4))

        # Number of elements, of nodes, of nodes per element and the magic number
        self.file.write(header.pack_int(4 * 4))
        self.file.write(header.pack_int(header.nb_elements))
        self.file.write(header.pack_int(header.nb_nodes))
        self.file.write(header.pack_int(header.nb_nodes_per_elem))
        self.file.write(header.pack_int(1))  # magic number
        self.file.write(header.pack_int(4 * 4))

        # IKLE
        nb_ikle_values = header.nb_elements * header.nb_nodes_per_elem
        self.file.write(header.pack_int(4 * nb_ikle_values))
        self.file.write(header.pack_int(*header.ikle, nb=nb_ikle_values))
        self.file.write(header.pack_int(4 * nb_ikle_values))

        # IPOBO
        self.file.write(header.pack_int(4 * header.nb_nodes))
        self.file.write(header.pack_int(*header.ipobo, nb=header.nb_nodes))
        self.file.write(header.pack_int(4 * header.nb_nodes))

        # X coordinates
        self.file.write(header.pack_int(header.float_size * header.nb_nodes))
        self.file.write(header.pack_float(*header.x_stored, nb=header.nb_nodes))
        self.file.write(header.pack_int(header.float_size * header.nb_nodes))

        # Y coordinates
        self.file.write(header.pack_int(header.float_size * header.nb_nodes))
        self.file.write(header.pack_float(*header.y_stored, nb=header.nb_nodes))
        self.file.write(header.pack_int(header.float_size * header.nb_nodes))

    def write_entire_frame(self, time_to_write, values):
        """!
        @brief write all variables/nodes values
        @param header <SerafinHeader>: output header
        @param time_to_write <float>: output time (in seconds)
        @param values <numpy 2D-array>: values to write, of dimension (nb_var, nb_nodes)
        """
        self.file.write(self.header.pack_int(header.float_size))
        self.file.write(header.pack_float(time_to_write))
        self.file.write(header.pack_int(header.float_size))

        for i in range(header.nb_var):
            self.file.write(header.pack_int(header.float_size * header.nb_nodes))
            self.file.write(header.pack_float(*values[i, :], nb=header.nb_nodes))
            self.file.write(header.pack_int(header.float_size * header.nb_nodes))
