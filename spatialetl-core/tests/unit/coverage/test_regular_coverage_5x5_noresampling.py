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
x_axis = np.arange(0, x_size)
y_axis = np.arange(0, y_size)


def test_axis_x(caplog):
    caplog.set_level(logging.INFO)

    for thread in range(2, os.cpu_count()):
        logging.info(f"Testing with {thread} threads")
        reader = MemoryReader(
            x=x_axis,
            y=y_axis
        )
        coverage = Coverage(reader=reader, nb_thread=thread)
        actual_data = coverage.read_axis_x()
        np.testing.assert_array_equal(actual_data, x_axis)


def test_axis_y(caplog):
    caplog.set_level(logging.INFO)

    for thread in range(2, os.cpu_count()):
        logging.info(f"Testing with {thread} threads")
        reader = MemoryReader(
            x=x_axis,
            y=y_axis
        )
        coverage = Coverage(reader=reader, nb_thread=thread)
        actual_data = coverage.read_axis_y()
        np.testing.assert_array_equal(actual_data, y_axis)


def test_bathymetry(caplog):
    caplog.set_level(logging.INFO)

    for thread in range(2, os.cpu_count()):
        logging.info(f"Testing with {thread} threads")
        reader = MemoryReader(
            x=x_axis,
            y=y_axis,
            bathy=data
        )
        coverage = Coverage(reader=reader, nb_thread=thread)
        actual_data = coverage.read_variable_bathymetry()
        np.testing.assert_array_equal(actual_data, data)
