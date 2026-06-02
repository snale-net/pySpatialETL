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

y_size = 10
x_size = 10
data = np.zeros([y_size, x_size])
np.fill_diagonal(data, 2)
source_x_axis = np.asarray([[2.3781743, 2.37824526, 2.37831553, 2.37838509, 2.37845394, 2.37852208
                                , 2.37858951, 2.37865624, 2.37872226, 2.37878756],
                            [2.37850636, 2.37857763, 2.37864818, 2.37871803, 2.37878718, 2.37885563
                                , 2.37892336, 2.37899039, 2.37905671, 2.37912232],
                            [2.37883913, 2.37891069, 2.37898154, 2.37905169, 2.37912114, 2.37918988
                                , 2.37925792, 2.37932525, 2.37939187, 2.3794578],
                            [2.3791726, 2.37924445, 2.3793156, 2.37938605, 2.3794558, 2.37952485
                                , 2.37959318, 2.37966082, 2.37972775, 2.37979397],
                            [2.37950676, 2.37957891, 2.37965036, 2.37972111, 2.37979117, 2.37986051
                                , 2.37992915, 2.37999709, 2.38006433, 2.38013086],
                            [2.37984163, 2.37991408, 2.37998583, 2.38005688, 2.38012723, 2.38019688
                                , 2.38026583, 2.38033407, 2.38040161, 2.38046845],
                            [2.38017719, 2.38024994, 2.38032199, 2.38039334, 2.380464, 2.38053395
                                , 2.3806032, 2.38067176, 2.3807396, 2.38080675],
                            [2.38051346, 2.3805865, 2.38065886, 2.38073051, 2.38080146, 2.38087172
                                , 2.38094128, 2.38101014, 2.38107829, 2.38114575],
                            [2.38085042, 2.38092376, 2.38099641, 2.38106837, 2.38113963, 2.38121019
                                , 2.38128006, 2.38134922, 2.38141769, 2.38148545],
                            [2.38118806, 2.38126171, 2.38133467, 2.38140692, 2.38147849, 2.38154936
                                , 2.38161953, 2.38168901, 2.38175778, 2.38182586]])
source_y_axis = np.asarray([[6.45008085, 6.4497511, 6.44942105, 6.4490907, 6.44876006, 6.44842912
                                , 6.44809788, 6.44776635, 6.44743451, 6.44710237],
                            [6.45015187, 6.44982141, 6.44949066, 6.44915961, 6.44882827, 6.44849663
                                , 6.44816469, 6.44783245, 6.44749991, 6.44716707],
                            [6.45022317, 6.44989202, 6.44956057, 6.44922883, 6.44889678, 6.44856444
                                , 6.4482318, 6.44789885, 6.44756561, 6.44723206],
                            [6.45029478, 6.44996292, 6.44963078, 6.44929833, 6.44896559, 6.44863255
                                , 6.4482992, 6.44796556, 6.44763161, 6.44729736],
                            [6.45036667, 6.45003412, 6.44970128, 6.44936814, 6.4490347, 6.44870096
                                , 6.44836691, 6.44803257, 6.44769792, 6.44736296],
                            [6.45043886, 6.45010562, 6.44977208, 6.44943824, 6.4491041, 6.44876966
                                , 6.44843492, 6.44809988, 6.44776453, 6.44742888],
                            [6.45051135, 6.45017741, 6.44984318, 6.44950864, 6.44917381, 6.44883868
                                , 6.44850324, 6.4481675, 6.44783145, 6.44749509],
                            [6.45058413, 6.4502495, 6.44991457, 6.44957935, 6.44924382, 6.44890799
                                , 6.44857186, 6.44823542, 6.44789867, 6.44756162],
                            [6.45065721, 6.45032189, 6.44998627, 6.44965035, 6.44931413, 6.44897761
                                , 6.44864078, 6.44830365, 6.4479662, 6.44762845],
                            [6.45073058, 6.45039458, 6.45005827, 6.44972166, 6.44938474, 6.44904753
                                , 6.44871001, 6.44837218, 6.44803404, 6.4476956]])
target_x_axis = [2.3781743, 2.3786743, 2.3791743000000003,
                 2.3796743000000005, 2.3801743000000006, 2.380674300000001,
                 2.381174300000001, 2.381674300000001]
target_y_axis = [6.44710237, 6.447602369999999, 6.448102369999999,
                 6.448602369999999, 6.4491023699999985, 6.449602369999998,
                 6.450102369999998, 6.450602369999998]
target_data = np.zeros([y_size - 2, x_size - 2])
np.fill_diagonal(target_data, 2)


def test_axis_x(caplog):
    caplog.set_level(logging.INFO)

    for thread in range(2, os.cpu_count()):
        logging.info(f"Testing with {thread} threads")
        reader = MemoryReader(
            x_axis=source_x_axis,
            y_axis=source_y_axis
        )

        coverage = Coverage(reader=reader, resolution_x=0.0005, resolution_y=0.0005, nb_thread=thread)
        actual_data = coverage.read_axis_x()
        np.testing.assert_array_equal(actual_data, target_x_axis)


def test_axis_y(caplog):
    caplog.set_level(logging.INFO)

    for thread in range(2, os.cpu_count()):
        logging.info(f"Testing with {thread} threads")
        reader = MemoryReader(
            x_axis=source_x_axis,
            y_axis=source_y_axis
        )

        coverage = Coverage(reader=reader, resolution_x=0.0005, resolution_y=0.0005, nb_thread=thread)
        actual_data = coverage.read_axis_y()
        np.testing.assert_array_equal(actual_data, target_y_axis)


def test_bathymetry(caplog):
    caplog.set_level(logging.INFO)

    for thread in range(2, os.cpu_count()):
        logging.info(f"Testing with {thread} threads")
        reader = MemoryReader(
            x_axis=source_x_axis,
            y_axis=source_y_axis,
            bathymetry=data
        )

        coverage = Coverage(reader=reader, resolution_x=0.0005, resolution_y=0.0005, nb_thread=thread)
        actual_data = coverage.read_variable_bathymetry()
        np.testing.assert_array_equal(actual_data, target_data)
