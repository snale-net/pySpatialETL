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
from concurrent.futures import ProcessPoolExecutor
import concurrent.futures

import numpy as np
import pytest
import multiprocessing
from multiprocessing import Pool
from spatialetl.coverage.coverage import Coverage
from spatialetl.coverage.io.memory_reader import MemoryReader
from spatialetl.providers.common.gdal.coverage.tiff.default_writer import DefaultWriter
from spatialetl.utils.logger import logging

def test_tiff_default_reader():
    data =np.zeros([10, 10])
    np.fill_diagonal(data, 5)
    reader = MemoryReader(
        x_axis=[0, 1, 2, 3, 4, 5, 6, 7, 8, 9],
        y_axis=[0, 1, 2, 3, 4, 5, 6, 7, 8, 9],
        bathymetry=data
    )

    logging.setLevel(logging.DEBUG)

    coverage = Coverage(reader=reader)

    writer = DefaultWriter(coverage,"/tmp")
    writer.write_variable_bathymetry()

    #np.testing.assert_array_equal( coverage.read_variable_bathymetry(),data)

    # Version process
    # nb_process = 2
    # for i in range(nb_process):
    #     coverage = Coverage(reader=reader, size=nb_process, rank=i)
    #
    #     p = multiprocessing.Process()
    #     p.start()
    #
    #     if coverage.rank == 0:
    #         assert coverage.read_axis_x() == [0, 1, 2, 3, 4, 5, 6, 7, 8, 9]
    #         assert coverage.read_axis_y() == [0, 1, 2, 3, 4]
    #         np.testing.assert_array_equal(coverage.read_variable_bathymetry(),[[5.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0], [0.0, 5.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0], [0.0, 0.0, 5.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0], [0.0, 0.0, 0.0, 5.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0], [0.0, 0.0, 0.0, 0.0, 5.0, 0.0, 0.0, 0.0, 0.0, 0.0]])
    #     elif coverage.rank == 1:
    #         assert coverage.read_axis_x() == [0, 1, 2, 3, 4, 5, 6, 7, 8, 9]
    #         assert coverage.read_axis_y() == [5, 6, 7, 8, 9]
    #         np.testing.assert_array_equal(coverage.read_variable_bathymetry(),[[0.0, 0.0, 0.0, 0.0, 0.0, 5.0, 0.0, 0.0, 0.0, 0.0], [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 5.0, 0.0, 0.0, 0.0], [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 5.0, 0.0, 0.0], [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 5.0, 0.0], [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 5.0]])
    #

