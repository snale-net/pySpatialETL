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
y_axis, x_axis = np.mgrid[0:y_size, 0:x_size]

@pytest.mark.mpi(ranks=[2,3], timeout=10, unit="s")
def test_axis_x(caplog, mpi_ranks):
    caplog.set_level(logging.INFO)

    for thread in range(2, os.cpu_count() - mpi_ranks):
        logging.info(f"Testing with {thread} threads")
        thread=1
        reader = MemoryReader(
            x=x_axis,
            y=y_axis
        )

        logging.setLevel(logging.DEBUG)

        coverage = Coverage(reader=reader, nb_thread=thread)

        if coverage.rank == 0:
            global_data = np.empty([coverage.get_y_size(type="target_global"),coverage.get_x_size(type="target_global")])
            global_data[:] = np.nan

        local_data = coverage.read_axis_x()

        if coverage.rank != 0:
            coverage.comm.Send(np.ascontiguousarray(local_data), dest=0)
        else:
            # Pour le proc n°1
            global_data[coverage.parallel_map[coverage.rank]["dst_global_y"],coverage.parallel_map[coverage.rank]["dst_global_x"]] = local_data

            # Pour les autres
            for source in range(1, coverage.size):
                recvbuf = np.empty([coverage.parallel_map[source]["dst_local_y_size"],coverage.parallel_map[source]["dst_local_x_size"]])
                coverage.comm.Recv(recvbuf, source=source)
                global_data[coverage.parallel_map[source]["dst_global_y"],coverage.parallel_map[source]["dst_global_x"]] = recvbuf

        coverage.comm.barrier()

        if coverage.rank == 0:
            np.testing.assert_array_equal(global_data, x_axis)

@pytest.mark.mpi(ranks=[2,3], timeout=10, unit="s")
def test_axis_y(caplog, mpi_ranks):
    caplog.set_level(logging.DEBUG)

    for thread in range(2, os.cpu_count() - mpi_ranks):
        logging.info(f"Testing with {thread} threads")
        reader = MemoryReader(
            x=x_axis,
            y=y_axis
        )

        coverage = Coverage(reader=reader,nb_thread=thread)

        if coverage.rank == 0:
            global_data = np.empty([coverage.get_y_size(type="target_global"),coverage.get_y_size(type="target_global")])
            global_data[:] = np.nan

        local_data = coverage.read_axis_y()

        if coverage.rank != 0:
            coverage.comm.Send(np.ascontiguousarray(local_data), dest=0)
        else:
            # Pour le proc n°1
            global_data[coverage.parallel_map[coverage.rank]["dst_global_y"],coverage.parallel_map[coverage.rank]["dst_global_y"]] = local_data

            # Pour les autres
            for source in range(1, coverage.size):
                recvbuf = np.empty([coverage.parallel_map[source]["dst_local_y_size"],coverage.parallel_map[source]["dst_local_x_size"]])
                coverage.comm.Recv(recvbuf, source=source)
                global_data[coverage.parallel_map[source]["dst_global_y"],coverage.parallel_map[source]["dst_global_x"]] = recvbuf

        coverage.comm.barrier()

        if coverage.rank == 0:
            np.testing.assert_array_equal(global_data,y_axis)


@pytest.mark.mpi(ranks=[2,3], timeout=10, unit="s")
def test_bathymetry(caplog,mpi_ranks):
    caplog.set_level(logging.DEBUG)

    for thread in range(2, os.cpu_count() - mpi_ranks):
        logging.info(f"Testing with {thread} threads")
        reader = MemoryReader(
            x=x_axis,
            y=y_axis,
            bathy=data
        )

        coverage = Coverage(reader=reader, nb_thread=thread)

        if coverage.rank == 0:
            global_data = np.empty(
                [coverage.get_y_size(type="target_global"),
                 coverage.get_x_size(type="target_global")])
            global_data[:] = np.nan

        local_data = coverage.read_variable_bathymetry()

        if coverage.rank != 0:
            coverage.comm.Send(np.ascontiguousarray(local_data), dest=0)
        else:
            # Pour le proc n°1
            global_data[coverage.parallel_map[coverage.rank]["dst_global_y"],
            coverage.parallel_map[coverage.rank]["dst_global_x"]] = local_data

            # Pour les autres
            for source in range(1, coverage.size):
                recvbuf = np.empty([coverage.parallel_map[source]["dst_local_y_size"],
                                    coverage.parallel_map[source]["dst_local_x_size"]])
                coverage.comm.Recv(recvbuf, source=source)

                global_data[coverage.parallel_map[source]["dst_global_y"],
                coverage.parallel_map[source]["dst_global_x"]] = recvbuf

        coverage.comm.barrier()

        if coverage.rank == 0:
            np.testing.assert_array_equal(global_data, data)
