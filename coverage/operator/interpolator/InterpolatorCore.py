#! /usr/bin/env python2.7
# -*- coding: utf-8 -*-
#
# CoverageProcessing is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# any later version.
#
# CoverageProcessing is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# Author : Fabien Rétif - fabien.retif@zoho.com
#
from __future__ import division, print_function, absolute_import
from scipy.interpolate import griddata
from scipy.interpolate import interp1d
import logging
import numpy as np

def resample_2d_to_grid(gridX,gridY,newX,newY,data,method):

    logging.debug("[InterpolatorCore][horizontal_interpolation()] starting interpolation with method '" + str(method) + "'")

    points = np.array([gridX.flatten(), gridY.flatten()]).T
    values = data.flatten()
    xx, yy = np.meshgrid(newX, newY)
    return griddata(points, values, (xx, yy), method=method, rescale=True,fill_value=9.96921e+36)

def vertical_interpolation(sourceAxis,targetAxis,data,method):
    logging.debug("[InterpolatorCore][vertical_interpolation()] Looking for water depth : " + str(
       targetAxis[0]) + " m with method '" + str(method) + "'.")
    #logging.debug("[InterpolatorCore][vertical_interpolation()] Source Axis contains: " + str(sourceAxis))
    #logging.debug("[InterpolatorCore][vertical_interpolation()] Candidates values are: " + str(data))
    #logging.debug("[InterpolatorCore][vertical_interpolation()] Target Axis contains: " + str(targetAxis))

    if method is None:
        return np.nan

    #if not np.isnan(data.all()):
    #    return np.nan

    elif method == "mean":
        return np.mean(data)

    elif method == "nearest":
        array = np.asarray(sourceAxis)
        nearest_index_t = (np.abs(array - targetAxis[0])).argmin()
        return data[nearest_index_t]

    else:
        f = interp1d(sourceAxis, data, kind=method, bounds_error=False)
        return f(targetAxis)
