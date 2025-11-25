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
import datetime

import numpy as np
import os
import pandas

from spatialetl.coverage import TimeCoverage
from spatialetl.coverage.coverage import Coverage
from spatialetl.coverage.io.memory_reader import MemoryReader
from spatialetl.utils.logger import logging

t_size = 10
y_size = 10
x_size = 10
data = np.ones([t_size,y_size, x_size])
#np.fill_diagonal(data, 2)
source_t_axis = pandas.date_range(start=datetime.datetime(2025,1,1,0,0,0),end=datetime.datetime(2025,1,1,9,0,0), freq="1H").to_pydatetime()
source_y_axis, source_x_axis = np.mgrid[0:y_size, 0:x_size]
target_x_axis = [0., 1., 2., 3., 4., 5., 6., 7., 8.]
target_y_axis = [0., 1., 2., 3., 4., 5., 6., 7., 8.]
target_data = np.ones([t_size,y_size - 1, x_size - 1])
#np.fill_diagonal(target_data, 2)

def test_axis_x(caplog):
    caplog.set_level(logging.INFO)

    for thread in range(2, os.cpu_count()):
        logging.info(f"Testing with {thread} threads")
        reader = MemoryReader(
            t=source_t_axis,
            x=source_x_axis,
            y=source_y_axis
        )

        coverage = TimeCoverage(reader=reader, resolution_x=1, resolution_y=1, nb_thread=thread)
        actual_data = coverage.read_axis_x()
        np.testing.assert_array_equal(actual_data, target_x_axis)


def test_axis_y(caplog):
    caplog.set_level(logging.INFO)

    for thread in range(2, os.cpu_count()):
        logging.info(f"Testing with {thread} threads")
        reader = MemoryReader(
            t=source_t_axis,
            x=source_x_axis,
            y=source_y_axis
        )

        coverage = TimeCoverage(reader=reader, resolution_x=1, resolution_y=1, nb_thread=thread)
        actual_data = coverage.read_axis_y()
        np.testing.assert_array_equal(actual_data, target_y_axis)


def test_surface_air_pressure(caplog):
    caplog.set_level(logging.DEBUG)

    for thread in range(2, os.cpu_count()):
        logging.info(f"Testing with {thread} threads")
        reader = MemoryReader(
            t=source_t_axis,
            x=source_x_axis,
            y=source_y_axis,
            sp=data
        )

        coverage = TimeCoverage(reader=reader, resolution_x=1, resolution_y=1, nb_thread=thread)
        index_time = 0
        for time in coverage.read_axis_t():
            actual_data = coverage.read_variable_surface_air_pressure_at_time(time)
            np.testing.assert_array_equal(actual_data, target_data[index_time,])
            index_time +=1
