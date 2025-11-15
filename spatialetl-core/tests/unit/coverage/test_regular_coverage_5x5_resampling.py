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

from spatialetl.coverage.coverage import Coverage
from spatialetl.coverage.io.memory_reader import MemoryReader
from spatialetl.utils.logger import logging

y_size = 5
x_size = 5
data = np.zeros([y_size, x_size])
np.fill_diagonal(data, 2)
source_x_axis = np.arange(0, x_size)
source_y_axis = np.arange(0, y_size)
target_x_axis = [0., 1., 2., 3.]
target_y_axis = [0., 1., 2., 3.]
target_data = np.zeros([y_size - 1, x_size - 1])
np.fill_diagonal(target_data, 2)


def test_axis_x(caplog):
    caplog.set_level(logging.INFO)

    for thread in range(2, os.cpu_count()):
        reader = MemoryReader(
            x=source_x_axis,
            y=source_y_axis
        )

        coverage = Coverage(reader=reader, resolution_x=1, resolution_y=1, nb_thread=thread)
        actual_data = coverage.read_axis_x()
        np.testing.assert_array_equal(actual_data, target_x_axis)


def test_axis_y(caplog):
    caplog.set_level(logging.INFO)

    for thread in range(2, os.cpu_count()):
        reader = MemoryReader(
            x=source_x_axis,
            y=source_y_axis
        )

        coverage = Coverage(reader=reader, resolution_x=1, resolution_y=1, nb_thread=thread)
        actual_data = coverage.read_axis_y()
        np.testing.assert_array_equal(actual_data, target_y_axis)


def test_bathymetry(caplog):
    caplog.set_level(logging.INFO)

    for thread in range(2, os.cpu_count()):
        reader = MemoryReader(
            x=source_x_axis,
            y=source_y_axis,
            bathy=data
        )

        coverage = Coverage(reader=reader, resolution_x=1, resolution_y=1, nb_thread=thread)
        actual_data = coverage.read_variable_bathymetry()
        np.testing.assert_array_equal(actual_data, target_data)
