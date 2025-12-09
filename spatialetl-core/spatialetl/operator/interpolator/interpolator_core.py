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
from __future__ import division, print_function, absolute_import

from scipy.spatial import qhull

from datetime import datetime

import numpy as np
from numpy import int8, int16, int32, int64
from scipy.interpolate import griddata, LinearNDInterpolator, NearestNDInterpolator
from scipy.interpolate import interp1d

from spatialetl.operator.parinterp.interpolator import Linear2DInterpolator
from spatialetl.utils.logger import logging


def interp_2d_weights(gridX, gridY, current_thread:int):
    """
    Interpolate the 2d weights of the source grid

    :param gridX:
    :param gridY:
    :param current_thread (int): Current thread
    :return current_thread, Delaunay triangulation
    """

    logging.debug(f"[InterpolatorCore][horizontal_interpolation()] Thread {current_thread} is computing weights")

    gridX = np.ma.filled(gridX, fill_value=-9999.)
    gridY = np.ma.filled(gridY, fill_value=-9999.)

    if gridX.ndim == 1 and gridY.ndim == 1:
        gridX, gridY = np.meshgrid(gridX, gridY)

    points = np.array([gridX.flatten(), gridY.flatten()]).T
    return current_thread, qhull.Delaunay(points)

def resample_2d_to_grid(gridX,gridY,newX,newY,data,method,current_thread):
    """
    2D resampling function
    """

    logging.debug(f"[InterpolatorCore][horizontal_interpolation()] Thread {current_thread} is starting interpolation with method '{method}'")

    gridX = np.ma.filled(gridX,fill_value=-9999.)
    gridY = np.ma.filled(gridY,fill_value=-9999.)

    if gridX.ndim ==1 and gridY.ndim==1:
        gridX, gridY = np.meshgrid(gridX, gridY)

    points = np.array([gridX.flatten(), gridY.flatten()]).T
    values = data.flatten()
    xx, yy = np.meshgrid(newX, newY)
    if data.dtype == int8 or data.dtype == int16 or data.dtype == int32 or data.dtype == int64:
        fill_value = -9999
    else:
        fill_value = 9.96921e+36

    #n, m = 300, 200
    #points = np.random.randint(low=0, high=1024, size=(n, m))
    #points = np.unique(points, axis=0)
    #x_points = points[: n // 2]

    #print(np.shape(x_points))
    #values = np.random.uniform(low=0.0, high=1.0, size=(len(x_points),))
    interp_points = np.array([xx.flatten(), yy.flatten()]).T

    print(np.shape(interp_points))

    interpet = Linear2DInterpolator(points, values)
    val = interpet(interp_points, values, fill_value=0.0)

    print(f"shape {np.shape(val)} = {val}")
    print(f"shape {np.shape(val.reshape(np.shape(newX)[0],np.shape(newY)[0]))} = {val.reshape(np.shape(newX)[0],np.shape(newY)[0])}")
    #val = griddata(points, values, (xx, yy), method=method, rescale=False,fill_value=fill_value)

    return current_thread,val.reshape(np.shape(newX)[0],np.shape(newY)[0])

def resample_faster_2d_to_grid(tri,newX,newY,data,method,current_thread):
    """
    2D resampling function
    """

    logging.debug(f"[InterpolatorCore][horizontal_interpolation()] Thread {current_thread} is starting interpolation")

    values = data.flatten()
    xx, yy = np.meshgrid(newX, newY)
    if data.dtype == int8 or data.dtype == int16 or data.dtype == int32 or data.dtype == int64:
        fill_value = -9999
    else:
        fill_value = 9.96921e+36

    if method =="linear":
        ip = LinearNDInterpolator(tri, values, fill_value=fill_value,
                              rescale=False)
    elif method == "nearest":
        ip = NearestNDInterpolator(tri, values, fill_value=fill_value,
                                  rescale=False)

    return current_thread,ip((xx, yy))


def vertical_interpolation(sourceAxis,targetAxis,data,method,extrapolate=False):
    #logging.debug("[InterpolatorCore][vertical_interpolation()] Looking for water depth : " + str(
    #   targetAxis[0]) + " m with method '" + str(method) + "'.")
    logging.debug("[InterpolatorCore][vertical_interpolation()] Source Axis contains: " + str(sourceAxis))
    logging.debug("[InterpolatorCore][vertical_interpolation()] Candidates values are: " + str(data))
    logging.debug("[InterpolatorCore][vertical_interpolation()] Target Axis contains: " + str(targetAxis))
    logging.debug("[InterpolatorCore][vertical_interpolation()] Method: " + str(method))
    logging.debug("[InterpolatorCore][vertical_interpolation()] ----------------------------------------")

    if method == "mean":
        return np.mean(data)

    elif method == "nearest":
        array = np.asarray(sourceAxis)
        nearest_index_t = (np.abs(array - targetAxis[0])).argmin()
        return data[nearest_index_t]
    elif method == "linear" or method == "cubic":
        try:
            f = interp1d(sourceAxis, data, kind=method, bounds_error=True)
            return f(targetAxis)
        except ValueError as ex:
            logging.warning("[InterpolatorCore][vertical_interpolation()] Error: " + str(ex))
            logging.warning("[InterpolatorCore][vertical_interpolation()] This error may occur when you ask for a water depth out of range: "+ str(targetAxis[0])+" m")
            logging.warning("[InterpolatorCore][vertical_interpolation()] We found these water depth candidates: " + str(sourceAxis))
            logging.warning("[InterpolatorCore][vertical_interpolation()] To avoid this error, you can change your zbox range or use another vertical interpolation method.")
            if extrapolate:
                logging.warning("[InterpolatorCore][vertical_interpolation()] We continue by using an extrapolation method.")
                logging.warning("[InterpolatorCore][vertical_interpolation()] ----------------------------------------")
                f = interp1d(sourceAxis, data, kind=method,  fill_value = "extrapolate")
                return f(targetAxis)
            else:
                logging.warning("[InterpolatorCore][vertical_interpolation()] We continue by using the nearest method.")
                logging.warning("[InterpolatorCore][vertical_interpolation()] ----------------------------------------")
                array = np.asarray(sourceAxis)
                nearest_index_t = (np.abs(array - targetAxis[0])).argmin()
                return data[nearest_index_t]
    else:
        raise ValueError("Unable to decode vertical interpolation method : "+str(method))

def time_1d_interpolation(sourceAxis,targetAxis,data,method,extrapolate=False):
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

    elif method == "linear" or method == "cubic":

        try:
            f = interp1d(sourceAxis, data, kind=method, bounds_error=True)
            return f(targetAxis)
        except ValueError as ex:
            logging.warning("[InterpolatorCore][time_interpolation()] Error: " + str(ex))
            logging.warning("[InterpolatorCore][time_interpolation()] This error may occur when you ask for a datetime  out of range : " + str(
                    datetime.utcfromtimestamp(targetAxis[0])))
            logging.warning("[InterpolatorCore][time_interpolation()] To avoid this error, you can change your time range or use another time interpolation method.")
            if extrapolate:
                logging.warning(
                    "[InterpolatorCore][vertical_interpolation()] We continue by using an extrapolation method.")
                logging.warning("[InterpolatorCore][vertical_interpolation()] ----------------------------------------")
                f = interp1d(sourceAxis, data, kind=method, fill_value="extrapolate")
                return f(targetAxis)
            else:
                logging.warning("[InterpolatorCore][vertical_interpolation()] We continue by using the nearest method.")
                logging.warning("[InterpolatorCore][vertical_interpolation()] ----------------------------------------")
                array = np.asarray(sourceAxis)
                nearest_index_t = (np.abs(array - targetAxis[0])).argmin()
                return data[nearest_index_t]
    else:
        raise ValueError("Unable to decode vertical interpolation method : " + str(method))

