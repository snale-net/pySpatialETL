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
import cftime
import numpy as np
import pytest
from spatialetl.providers.symphonie.coverage.netcdf.symphonie.v293.symphonie_reader import SYMPHONIEReader

from spatialetl.coverage import TimeCoverage
from spatialetl.exception.not_found_in_rank_error import NotFoundInRankError

def test_mpi_coverage():
    reader = SYMPHONIEReader("../io/netcdf/symphonie/v293/tests/resources/grid.nc",
                             "../io/netcdf/symphonie/v293/tests/resources/2014*")
    coverage = TimeCoverage(reader)

    # test_get_t_size()
    expected_size = coverage.parallel_map[coverage.rank]["dst_local_t_size"]
    candidate_size = coverage.get_t_size()
    assert expected_size==candidate_size, "test_get_t_size()"

    # test_read_axis_t()
    expected_shape = (coverage.parallel_map[coverage.rank]["dst_local_t_size"])
    candidate_shape = np.shape(coverage.read_axis_t())
    assert expected_shape==candidate_shape, "test_read_axis_t()"

    # test_find_time_index()
    with pytest.raises(NotFoundInRankError):
        expected_value = 2
        date = cftime.datetime(2014, 1, 1, 13, 0, 00)
        candidate_value = coverage.find_time_index(date, domain="source_global")
        assert expected_value==candidate_value, "test_find_time_index()"


    # test_read_variable_bathymetry()
    # expected_shape = (coverage.parallel_map[coverage.rank]["dst_local_y_size"],
    #                  coverage.parallel_map[coverage.rank]["dst_local_x_size"])
    # candidate_shape = np.shape(coverage.read_variable_bathymetry())
    # assert expected_shape, candidate_shape, "test_read_variable_bathymetry()")

def test_mpi_interpolated_coverage():
    reader = SYMPHONIEReader("../io/netcdf/symphonie/v293/tests/resources/grid.nc",
                             "../io/netcdf/symphonie/v293/tests/resources/2014*")
    coverage = TimeCoverage(reader, freq="45min")

    # test_get_t_size()
    expected_size = coverage.parallel_map[coverage.rank]["dst_local_t_size"]
    candidate_size = coverage.get_t_size()
    assert expected_size==candidate_size, "test_get_t_size()"

    # test_read_axis_t()
    expected_shape = (coverage.parallel_map[coverage.rank]["dst_local_t_size"])
    candidate_shape = np.shape(coverage.read_axis_t())
    assert expected_shape==candidate_shape, "test_read_axis_t()"

    # test_find_time_index()
    with pytest.raises(NotFoundInRankError):
        expected_value = 2
        date = cftime.datetime(2014, 1, 1, 13, 0, 00)
        # print(coverage.read_axis_t(type="source",timestamp=1))
        candidate_value = coverage.find_time_index(date, domain="source_global")
        assert expected_value==candidate_value, "test_find_time_index()"

