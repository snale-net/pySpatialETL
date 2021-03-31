#! /usr/bin/env python2.7
# -*- coding: utf-8 -*-
#
# pySpatialETL is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# any later version.
#
# pySpatialETL is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#

# Lien vers le dossier de la lib
import sys
sys.path = ['/work/sciences/toolbox/python/pyGeoSpatialETL'] + sys.path

from point.TimeMultiPoint import TimeMultiPoint
from point.io.netcdf.ecmwf.ECMWFReader import ECMWFReader
import logging
from datetime import timedelta
from point.io.netcdf.DefaultTimeMultiPointWriter import DefaultTimeMultiPointWriter as NcWriter
from point.io.ascii.DefaultTimePointWriter import DefaultTimePointWriter as AsciiWriter


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

    print('End of programm')

