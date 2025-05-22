#! /usr/bin/env python2.7
# -*- coding: utf-8 -*-
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
from spatialetl.point.TimeMultiPoint import TimeMultiPoint
from spatialetl.point.io.netcdf.ecmwf.ECMWFReader import ECMWFReader
from spatialetl.point.io.netcdf.DefaultWriter import DefaultWriter as NcWriter
from spatialetl.point.io.ascii.default_time_point_writer import DefaultTimePointWriter as AsciiWriter
import logging
from datetime import timedelta

if __name__ == "__main__":

    logging.basicConfig(format='[%(levelname)s] %(message)s', level=logging.INFO)

    points = {};
    points['IRHOB'] =   [2.43931,  6.36887]
    points['NOKOUE'] = [2.378550,   6.445640]
    stationCoords = [points['IRHOB'], points['NOKOUE']]

    TimeMultiPoint.TIME_DELTA_MIN = timedelta(hours=3)
    TimeMultiPoint.TIME_DELTA_MAX = timedelta(hours=3)
    TimeMultiPoint.TIME_INTERPOLATION_METHOD = "linear"

    reader = ECMWFReader('/work/sciences/projects/Nokoue-2019/forcings/ecmwf/previ-Benin.2018.nc',
                         stationCoords)

    myPoints = TimeMultiPoint(reader,start='2018-01-02T00:00:00',
                              end='2018-02-01T00:00:00', freq='3H')

    # NetCDF export
    writer = NcWriter(myPoints, '/tmp/ecmwf_points.nc')
    writer.write_variable_wind_speed_10m()
    writer.write_variable_wind_from_direction_10m()
    writer.write_variable_wind_to_direction_10m()
    writer.write_variable_surface_air_pressure()
    writer.write_variable_rainfall_amount()
    writer.close()

    # ASCII export
    for i in range(0,myPoints.get_nb_points()):
        writer = AsciiWriter(myPoints,i, '/tmp/ecmwf_points_'+str(i)+'.dat')
        writer.write_variable_wind_speed_10m()
        writer.write_variable_wind_from_direction_10m()
        writer.write_variable_surface_air_pressure()
        writer.write_variable_rainfall_amount()
        writer.close()

    print('End of program')

