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
from numpy import int8,int16,int32,int64
from datetime import datetime

def resample_2d_to_grid(gridX,gridY,newX,newY,data,method):

    logging.debug("[InterpolatorCore][horizontal_interpolation()] starting interpolation with method '" + str(method) + "'")

    if gridX.ndim ==1 and gridY.ndim==1:
        gridX, gridY = np.meshgrid(gridX, gridY)

    points = np.array([gridX.flatten(), gridY.flatten()]).T
    values = data.flatten()
    xx, yy = np.meshgrid(newX, newY)
    if data.dtype == int8 or data.dtype == int16 or data.dtype == int32 or data.dtype == int64:
        fill_value = -9999
    else:
        fill_value = 9.96921e+36

    return griddata(points, values, (xx, yy), method=method, rescale=True,fill_value=fill_value)

def vertical_interpolation(sourceAxis,targetAxis,data,method):
    #logging.debug("[InterpolatorCore][vertical_interpolation()] Looking for water depth : " + str(
    #   targetAxis[0]) + " m with method '" + str(method) + "'.")
    logging.debug("[InterpolatorCore][vertical_interpolation()] Source Axis contains: " + str(sourceAxis))
    logging.debug("[InterpolatorCore][vertical_interpolation()] Candidates values are: " + str(data))
    logging.debug("[InterpolatorCore][vertical_interpolation()] Target Axis contains: " + str(targetAxis))
    logging.debug("[InterpolatorCore][vertical_interpolation()] Method: " + str(method))
    logging.debug("[InterpolatorCore][vertical_interpolation()] ----------------------------------------")

    if method is None:
        return np.nan

    elif method == "mean":
        return np.mean(data)

    elif method == "nearest":
        array = np.asarray(sourceAxis)
        nearest_index_t = (np.abs(array - targetAxis[0])).argmin()
        return data[nearest_index_t]

    else:
        try:
            f = interp1d(sourceAxis, data, kind=method, bounds_error=True)
            return f(targetAxis)
        except ValueError as ex:
            logging.warning("[InterpolatorCore][vertical_interpolation()] Error: " + str(ex))
            logging.warning("[InterpolatorCore][vertical_interpolation()] This error may occur when you ask for a water depth out of range : "+ str(targetAxis[0])+" m")
            logging.warning("[InterpolatorCore][vertical_interpolation()] We continue by using an extrapolation method.")
            logging.warning("[InterpolatorCore][vertical_interpolation()] To avoid this error, you can change your zbox range or use another vertical interpolation method.")
            logging.warning("[InterpolatorCore][vertical_interpolation()] ----------------------------------------")
            f = interp1d(sourceAxis, data, kind=method,  fill_value = "extrapolate")
            return f(targetAxis)

def time_1d_interpolation(sourceAxis,targetAxis,data,method):
    logging.debug("[InterpolatorCore][time_interpolation()] Looking for time : "+str(datetime.utcfromtimestamp(targetAxis[0]))+" with method '"+str(method)+"'.")
    for time in sourceAxis:
        logging.debug("[InterpolatorCore][time_interpolation()] Source Axis contains: "+str(datetime.utcfromtimestamp(time)))

    for time in targetAxis:
        logging.debug("[InterpolatorCore][time_interpolation()] Target Axis contains: "+str(datetime.utcfromtimestamp(time)))

    if method is None:
        return np.nan

    elif method == "mean":
        return np.mean(data)

    elif method == "nearest":

        array = np.asarray(sourceAxis)
        nearest_index_t = (np.abs(array - targetAxis[0])).argmin()
        return data[nearest_index_t]

    else:
        f = interp1d(sourceAxis, data, kind=method, bounds_error=True)
        return f(targetAxis)

