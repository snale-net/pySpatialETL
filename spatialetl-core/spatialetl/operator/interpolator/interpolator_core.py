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

from datetime import datetime

import numpy as np
from numpy import int8, int16, int32, int64
from omp4py import omp, omp_get_thread_num
from scipy.interpolate import griddata, LinearNDInterpolator, NearestNDInterpolator
from scipy.interpolate import interp1d
from scipy.spatial import qhull

from spatialetl.operator.parinterp.interpolator import Linear2DInterpolator
from spatialetl.utils.logger import logging


# try:
#    from spatialetl.operator.parinterp.interpolator import Linear2DInterpolator
# except Exception as ex:
#    print(ex)

@omp()
def interp_2d_weights(gridX, gridY, map,threads_count:int):
    """
    Interpolate the 2d weights of the source grid

    :param gridX:
    :param gridY:
    :param current_thread (int): Current thread
    :return current_thread, Delaunay triangulation
    """

    logging.debug(f"[InterpolatorCore][horizontal_interpolation()] computing weights")

    gridX = np.ma.filled(gridX, fill_value=-9999.)
    gridY = np.ma.filled(gridY, fill_value=-9999.)

    tri = np.empty([threads_count], dtype=object)

    with omp("parallel shared(tri)"):
        current_thread = omp_get_thread_num()

        if gridX.ndim == 1 and gridY.ndim == 1:
            local_grid_x,  local_grid_y = np.meshgrid(gridX[map[current_thread]["src_global_x_overlap"]], gridY[map[current_thread]["src_global_y_overlap"]])
        else:
            local_grid_x = gridX[map[current_thread]["src_global_y_overlap"],
            map[current_thread]["src_global_x_overlap"]]
            local_grid_y = gridY[map[current_thread]["src_global_y_overlap"],
            map[current_thread]["src_global_x_overlap"]]

        points = np.array([local_grid_x.flatten(), local_grid_y.flatten()]).T
        tri[current_thread] = qhull.Delaunay(points)

    return tri

@omp()
def resample_2d_to_grid(gridX,gridY,newX,newY,data,map,method):
    """
    2D resampling function
    """

    logging.debug(f"[InterpolatorCore][horizontal_interpolation()] Thread {omp_get_thread_num()} is starting interpolation with method '{method}'")

    local_data = np.zeros([np.shape(newY)[0],np.shape(newX)[0]])

    with omp("parallel shared(local_data)"):
        current_thread = omp_get_thread_num()

        # source grid
        if gridX.ndim == 1 and gridY.ndim == 1:
            local_grid_x, local_grid_y = np.meshgrid(gridX[map[current_thread]["src_global_x_overlap"]],
                                                     gridY[map[current_thread]["src_global_y_overlap"]])
        else:
            local_grid_x = gridX[map[current_thread]["src_global_y_overlap"],
            map[current_thread]["src_global_x_overlap"]]
            local_grid_y = gridY[map[current_thread]["src_global_y_overlap"],
            map[current_thread]["src_global_x_overlap"]]

        points = np.array([local_grid_x.flatten(), local_grid_y.flatten()]).T

        # Target grid
        dst_local_grid_x, dst_local_grid_y = np.meshgrid(newX[map[current_thread]["dst_global_x_overlap"]],
                                                         newY[map[current_thread]["dst_global_y_overlap"]])

        # Values
        values = data[map[current_thread]["src_global_y_overlap"],
            map[current_thread]["src_global_x_overlap"]].flatten()


        if data.dtype == int8 or data.dtype == int16 or data.dtype == int32 or data.dtype == int64:
            fill_value = -9999
        else:
            fill_value = 9.96921e+36

        # scipy
        #data = griddata(points, values, (dst_local_grid_x, dst_local_grid_y), method=method, rescale=False, fill_value=fill_value)

        # parinterp
        interp_points = np.array([dst_local_grid_x.flatten(), dst_local_grid_y.flatten()]).T
        interpet = Linear2DInterpolator(points, values)
        data = interpet(interp_points, values, fill_value=0.0).reshape(np.shape(newX)[0],np.shape(newY)[0])

        local_data[map[current_thread]["dst_mpi_y"], map[current_thread]["dst_mpi_x"]] = data[
            map[current_thread]["dst_local_y"],
            map[current_thread]["dst_local_x"]]

    return local_data

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

