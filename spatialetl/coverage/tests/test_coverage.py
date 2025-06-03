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
import sys
sys.path = ['/work/sciences/pySpatialETL'] + sys.path
from unittest import TestCase
import numpy as np

from spatialetl.coverage import Coverage
from spatialetl.coverage.io.netcdf.symphonie.v293 import SYMPHONIEReader
from spatialetl.exception.not_found_in_rank_error import NotFoundInRankError


class TestCoverage(TestCase):

    def test_mpi_coverage(self):
        reader = SYMPHONIEReader("../io/netcdf/symphonie/v293/tests/resources/grid.nc", "../io/netcdf/symphonie/v293/tests/resources/2014*")
        coverage = Coverage(reader)

        # test_get_x_size()
        expected_size = coverage.map_mpi[coverage.rank]["dst_local_x_size"]
        candidate_size = coverage.get_x_size()
        self.assertEqual(expected_size, candidate_size, "test_get_x_size()")

        # test_get_y_size()
        expected_size = coverage.map_mpi[coverage.rank]["dst_local_y_size"]
        candidate_size = coverage.get_y_size()
        self.assertEqual(expected_size, candidate_size, "test_get_y_size()")

        # test_is_regular_grid()
        expected_value = False
        candidate_value = coverage.is_regular_grid()
        self.assertEqual(expected_value, candidate_value, "test_is_regular_grid()")

        # test_read_axis_x()
        expected_shape = (coverage.map_mpi[coverage.rank]["dst_local_y_size"],
                          coverage.map_mpi[coverage.rank]["dst_local_x_size"])
        candidate_shape = np.shape(coverage.read_axis_x())
        self.assertEqual(expected_shape, candidate_shape, "test_read_axis_x()")

        # test_read_axis_y()
        expected_shape = (coverage.map_mpi[coverage.rank]["dst_local_y_size"],
                          coverage.map_mpi[coverage.rank]["dst_local_x_size"])
        candidate_shape = np.shape(coverage.read_axis_y())
        self.assertEqual(expected_shape, candidate_shape, "test_read_axis_y()")

        # test_find_point_index()
        try:
            expected_value = [5, 6, 0.07194630440530268, 0.08993284357837819, 0.006314286576960388]
            candidate_value = coverage.find_point_index(0.0719, 0.0899, decimal_tolerance=3,type="source_global")
            self.assertEqual(expected_value, candidate_value, "test_find_point_index()")
        except NotFoundInRankError as ex:
            if coverage.size == 1:
                self.fail(ex.message)

        # test_read_variable_bathymetry()
        expected_shape = (coverage.map_mpi[coverage.rank]["dst_local_y_size"],
                          coverage.map_mpi[coverage.rank]["dst_local_x_size"])
        candidate_shape = np.shape(coverage.read_variable_bathymetry())
        self.assertEqual(expected_shape, candidate_shape, "test_read_variable_bathymetry()")

        # test_read_variable_mesh_size()
        expected_shape = (coverage.map_mpi[coverage.rank]["dst_local_y_size"],
                          coverage.map_mpi[coverage.rank]["dst_local_x_size"])
        candidate_shape = np.shape(coverage.read_variable_mesh_size())
        self.assertEqual(expected_shape, candidate_shape, "test_read_variable_mesh_size()")

        # test_read_variable_x_mesh_size()
        expected_shape = (coverage.map_mpi[coverage.rank]["dst_local_y_size"],
                          coverage.map_mpi[coverage.rank]["dst_local_x_size"])
        candidate_shape = np.shape(coverage.read_variable_x_mesh_size())
        self.assertEqual(expected_shape, candidate_shape, "test_read_variable_x_mesh_size()")

        # test_read_variable_y_mesh_size()
        expected_shape = (coverage.map_mpi[coverage.rank]["dst_local_y_size"],
                          coverage.map_mpi[coverage.rank]["dst_local_x_size"])
        candidate_shape = np.shape(coverage.read_variable_y_mesh_size())
        self.assertEqual(expected_shape, candidate_shape, "test_read_variable_y_mesh_size()")

        # test_read_variable_2D_sea_binary_mask()
        expected_shape = (coverage.map_mpi[coverage.rank]["dst_local_y_size"],
                          coverage.map_mpi[coverage.rank]["dst_local_x_size"])
        candidate_shape = np.shape(coverage.read_variable_2D_sea_binary_mask())
        self.assertEqual(expected_shape, candidate_shape, "test_read_variable_2D_sea_binary_mask()")

    def test_mpi_interpolated_coverage(self):
        reader = SYMPHONIEReader("../io/netcdf/symphonie/v293/tests/resources/grid.nc",
                                 "../io/netcdf/symphonie/v293/tests/resources/2014*")
        coverage = Coverage(reader, resolution_x=0.001, resolution_y=0.001)

        # test_get_x_size()
        expected_size = coverage.map_mpi[coverage.rank]["dst_local_x_size"]
        candidate_size = coverage.get_x_size()
        self.assertEqual(expected_size, candidate_size, "test_get_x_size()")

        # test_get_y_size()
        expected_size = coverage.map_mpi[coverage.rank]["dst_local_y_size"]
        candidate_size = coverage.get_y_size()
        self.assertEqual(expected_size, candidate_size, "test_get_y_size()")

        # test_is_regular_grid()
        expected_value = True
        candidate_value = coverage.is_regular_grid()
        self.assertEqual(expected_value, candidate_value, "test_is_regular_grid()")

        # test_read_axis_x()
        expected_shape = (coverage.map_mpi[coverage.rank]["dst_local_x_size"])
        candidate_shape = np.shape(coverage.read_axis_x())
        self.assertEqual(expected_shape, candidate_shape, "test_read_axis_x()")

        # test_read_axis_y()
        expected_shape = (coverage.map_mpi[coverage.rank]["dst_local_y_size"])
        candidate_shape = np.shape(coverage.read_axis_y())
        self.assertEqual(expected_shape, candidate_shape, "test_read_axis_y()")

        # test_find_point_index()
        try:
            expected_value = [5, 6, 0.07194630440530268, 0.08993284357837819, 0.006314286576960388]
            candidate_value = coverage.find_point_index(0.0719, 0.0899, decimal_tolerance=3, type="source_global")
            self.assertEqual(expected_value, candidate_value, "test_find_point_index()")
        except NotFoundInRankError as ex:
            if coverage.size == 1:
                self.fail(ex.message)

        # test_read_variable_bathymetry()
        expected_shape = (coverage.map_mpi[coverage.rank]["dst_local_y_size"],
                          coverage.map_mpi[coverage.rank]["dst_local_x_size"])
        candidate_shape = np.shape(coverage.read_variable_bathymetry())
        self.assertEqual(expected_shape, candidate_shape, "test_read_variable_bathymetry()")

        # test_read_variable_mesh_size()
        expected_shape = (coverage.map_mpi[coverage.rank]["dst_local_y_size"],
                          coverage.map_mpi[coverage.rank]["dst_local_x_size"])
        candidate_shape = np.shape(coverage.read_variable_mesh_size())
        self.assertEqual(expected_shape, candidate_shape, "test_read_variable_mesh_size()")

        # test_read_variable_x_mesh_size()
        expected_shape = (coverage.map_mpi[coverage.rank]["dst_local_y_size"],
                          coverage.map_mpi[coverage.rank]["dst_local_x_size"])
        candidate_shape = np.shape(coverage.read_variable_x_mesh_size())
        self.assertEqual(expected_shape, candidate_shape, "test_read_variable_x_mesh_size()")

        # test_read_variable_y_mesh_size()
        expected_shape = (coverage.map_mpi[coverage.rank]["dst_local_y_size"],
                          coverage.map_mpi[coverage.rank]["dst_local_x_size"])
        candidate_shape = np.shape(coverage.read_variable_y_mesh_size())
        self.assertEqual(expected_shape, candidate_shape, "test_read_variable_y_mesh_size()")

        # test_read_variable_2D_sea_binary_mask()
        expected_shape = (coverage.map_mpi[coverage.rank]["dst_local_y_size"],
                          coverage.map_mpi[coverage.rank]["dst_local_x_size"])
        candidate_shape = np.shape(coverage.read_variable_2D_sea_binary_mask())
        self.assertEqual(expected_shape, candidate_shape, "test_read_variable_2D_sea_binary_mask()")

