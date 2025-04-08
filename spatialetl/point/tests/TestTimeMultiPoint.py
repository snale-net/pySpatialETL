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
from unittest import TestCase
from datetime import datetime
from spatialetl.point.TimeMultiPoint import TimeMultiPoint
from spatialetl.point.io.ascii.DefaultTimePointReader import DefaultTimePointReader


class TestTimeMultiPoint(TestCase):

    def test_timeseries(self):
        reader = DefaultTimePointReader("../io/ascii/tests/ressources/Port-Sainte-Marie-la-Mer_03-2013_to_03-2013.dat",
                                        colsNumber=[0,1],
                                        varNames=["date","sea_water_column_thickness_at_time"])
        point = TimeMultiPoint(reader);

        # test_get_nb_points
        expected_size = 1
        candidate_size = point.get_nb_points()
        self.assertEqual(expected_size, candidate_size, "get_nb_points()")

        # test_read_variable_point_names
        expected_value = ["Port-Sainte-Marie-la-Mer"]
        candidate_value = point.read_variable_point_names()
        self.assertEqual(expected_value, candidate_value, "read_variable_point_names()")

        # test_read_x_axis
        expected_value = [3.0381070436335205]
        candidate_value = point.read_axis_x()
        self.assertEqual(expected_value, candidate_value, "read_axis_x()")

        # test_read_y_axis
        expected_value = [42.7230238237492]
        candidate_value = point.read_axis_y()
        self.assertEqual(expected_value, candidate_value, "read_axis_y()")

        # test_get_t_size
        expected_size = 13
        candidate_size = point.get_t_size()
        self.assertEqual(expected_size, candidate_size, "get_t_size()")

        # test_read_axis_t
        expected_value = datetime.strptime("2013-03-05 02:45:36", '%Y-%m-%d %H:%M:%S')
        candidate_value = point.read_axis_t()[9]
        self.assertEqual(expected_value, candidate_value, "read_axis_t()")

        # test_read_variable_sea_water_column_thickness_at_time
        expected_value = [1.2138863801956177]
        candidate_value = point.read_variable_sea_water_column_thickness_at_time(datetime.strptime("2013-03-05 02:45:36",'%Y-%m-%d %H:%M:%S'))
        self.assertEqual(expected_value, candidate_value, "read_variable_sea_water_column_thickness_at_time()")



