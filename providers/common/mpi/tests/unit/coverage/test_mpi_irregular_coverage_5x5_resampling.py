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
import numpy as np
import os

import pytest

from spatialetl.coverage.coverage import Coverage
from spatialetl.coverage.io.memory_reader import MemoryReader
from spatialetl.utils.logger import logging

y_size = 5
x_size = 5
data = np.zeros([y_size, x_size])
np.fill_diagonal(data, 2)
source_x_axis = np.asarray([[2.3781743, 2.37824526, 2.37831553, 2.37838509, 2.37845394],
                            [2.37850636, 2.37857763, 2.37864818, 2.37871803, 2.37878718],
                            [2.37883913, 2.37891069, 2.37898154, 2.37905169, 2.37912114],
                            [2.3791726, 2.37924445, 2.3793156, 2.37938605, 2.3794558],
                            [2.37950676, 2.37957891, 2.37965036, 2.37972111, 2.37979117]])
source_y_axis = np.asarray([[6.45008085, 6.4497511, 6.44942105, 6.4490907, 6.44876006],
                            [6.45015187, 6.44982141, 6.44949066, 6.44915961, 6.44882827],
                            [6.45022317, 6.44989202, 6.44956057, 6.44922883, 6.44889678],
                            [6.45029478, 6.44996292, 6.44963078, 6.44929833, 6.44896559],
                            [6.45036667, 6.45003412, 6.44970128, 6.44936814, 6.4490347]])
target_x_axis = [2.3781743, 2.3786743, 2.3791743000000003, 2.3796743000000005]
target_y_axis = [6.44876006, 6.449260059999999, 6.449760059999999, 6.450260059999999]
target_data = np.zeros([y_size - 2, x_size - 2])
np.fill_diagonal(target_data, 2)

@pytest.mark.mpi(ranks=[2,3], timeout=10, unit="s")
def test_axis_x(caplog, mpi_ranks):
    caplog.set_level(logging.INFO)

    for thread in range(2, os.cpu_count() - mpi_ranks):
        logging.info(f"Testing with {thread} threads")
        thread=1
        reader = MemoryReader(
            x_axis=source_x_axis,
            y_axis=source_y_axis
        )

        coverage = Coverage(reader=reader,resolution_x=0.0005, resolution_y=0.0005, nb_thread=thread)

        if coverage.rank == 0:
            actual_data = np.empty([coverage.get_y_size(type="target_global"),coverage.get_x_size(type="target_global")])
            actual_data[:] = np.nan

        local_data = coverage.read_axis_x()

        if coverage.rank != 0:
            coverage.comm.Send(np.ascontiguousarray(local_data), dest=0)
        else:
            # Pour le proc n°1
            actual_data[coverage.parallel_map[coverage.rank]["dst_global_y"],coverage.parallel_map[coverage.rank]["dst_global_x"]] = local_data

            # Pour les autres
            for source in range(1, coverage.size):
                recvbuf = np.empty([coverage.parallel_map[source]["dst_local_y_size"],coverage.parallel_map[source]["dst_local_x_size"]])
                coverage.comm.Recv(recvbuf, source=source)
                actual_data[coverage.parallel_map[source]["dst_global_y"],coverage.parallel_map[source]["dst_global_x"]] = recvbuf

        coverage.comm.barrier()

        if coverage.rank == 0:
            np.testing.assert_array_equal(actual_data, target_x_axis)

@pytest.mark.mpi(ranks=[2,3], timeout=10, unit="s")
def test_axis_y(caplog, mpi_ranks):
    caplog.set_level(logging.DEBUG)

    for thread in range(2, os.cpu_count() - mpi_ranks):
        logging.info(f"Testing with {thread} threads")
        reader = MemoryReader(
            x_axis=source_x_axis,
            y_axis=source_y_axis
        )

        coverage = Coverage(reader=reader,resolution_x=0.0005, resolution_y=0.0005,nb_thread=thread)

        if coverage.rank == 0:
            actual_data = np.empty([coverage.get_y_size(type="target_global"),coverage.get_y_size(type="target_global")])
            actual_data[:] = np.nan

        local_data = coverage.read_axis_y()

        if coverage.rank != 0:
            coverage.comm.Send(np.ascontiguousarray(local_data), dest=0)
        else:
            # Pour le proc n°1
            actual_data[coverage.parallel_map[coverage.rank]["dst_global_y"],coverage.parallel_map[coverage.rank]["dst_global_y"]] = local_data

            # Pour les autres
            for source in range(1, coverage.size):
                recvbuf = np.empty([coverage.parallel_map[source]["dst_local_y_size"],coverage.parallel_map[source]["dst_local_x_size"]])
                coverage.comm.Recv(recvbuf, source=source)
                actual_data[coverage.parallel_map[source]["dst_global_y"],coverage.parallel_map[source]["dst_global_x"]] = recvbuf

        coverage.comm.barrier()

        if coverage.rank == 0:
            np.testing.assert_array_equal(actual_data,target_y_axis)


@pytest.mark.mpi(ranks=[2,3], timeout=10, unit="s")
def test_bathymetry(caplog,mpi_ranks):
#def test_bathymetry():
    #caplog.set_level(logging.DEBUG)
    logging.setLevel(logging.DEBUG)

    for thread in range(2, os.cpu_count() - mpi_ranks):
        logging.info(f"Testing with {thread} threads")
        reader = MemoryReader(
            x_axis=source_x_axis,
            y_axis=source_y_axis,
            bathymetry=data
        )

        coverage = Coverage(reader=reader,resolution_x=0.0005, resolution_y=0.0005, nb_thread=thread)

        if coverage.rank == 0:
            actual_data = np.empty(
                [coverage.get_y_size(type="target_global"),
                 coverage.get_x_size(type="target_global")])
            actual_data[:] = np.nan

        local_data = coverage.read_variable_bathymetry()

        if coverage.rank != 0:
            coverage.comm.Send(np.ascontiguousarray(local_data), dest=0)
        else:
            # Pour le proc n°1
            actual_data[coverage.parallel_map[coverage.rank]["dst_global_y"],
            coverage.parallel_map[coverage.rank]["dst_global_x"]] = local_data

            # Pour les autres
            for source in range(1, coverage.size):
                recvbuf = np.empty([coverage.parallel_map[source]["dst_local_y_size"],
                                    coverage.parallel_map[source]["dst_local_x_size"]])
                coverage.comm.Recv(recvbuf, source=source)

                actual_data[coverage.parallel_map[source]["dst_global_y"],
                coverage.parallel_map[source]["dst_global_x"]] = recvbuf

        coverage.comm.barrier()

        if coverage.rank == 0:
            np.testing.assert_array_equal(actual_data, target_data)


if __name__ == '__main__':
    test_bathymetry()
