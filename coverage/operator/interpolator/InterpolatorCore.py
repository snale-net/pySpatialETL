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

class InterpolatorCore(object):

    HORIZONTAL_INTERPOLATION_METHOD = "nearest";
    #HORIZONTAL_INTERPOLATION_METHOD = "linear";
    VERTICAL_INTERPOLATION_METHOD = "nearest";

def resample_2d_to_grid(sourceAxisX,sourceAxisY,targetAxisX, targetAixsY,data):

    return griddata((sourceAxisX.ravel(),sourceAxisY.ravel()), data.ravel(),(targetAxisX, targetAixsY), method=InterpolatorCore.HORIZONTAL_INTERPOLATION_METHOD, fill_value=9.96921e+36)

def resample_2d_to_resolution(sourceAxisX,sourceAxisY,targetResX, targetResY,data):
    
    Ymin=np.min(sourceAxisY)
    Ymax=np.max(sourceAxisY)
    Xmin=np.min(sourceAxisX)
    Xmax=np.max(sourceAxisX)

    res=np.mean([targetResX,targetResY])
    
    lon_reg,lat_reg=np.meshgrid(np.arange(Xmin, Xmax, res),np.arange(Ymin, Ymax, res))

    return griddata((sourceAxisX.ravel(),sourceAxisY.ravel()), data.ravel(),(lon_reg, lat_reg), method=InterpolatorCore.HORIZONTAL_INTERPOLATION_METHOD)

def resample_2d(sourceAxisX,sourceAxisY,data):
    
    size=np.shape(data)    
    Ymin=np.min(sourceAxisY)
    Ymax=np.max(sourceAxisY)
    Xmin=np.min(sourceAxisX)
    Xmax=np.max(sourceAxisX)

    xres_mean=((Xmax-Xmin)/size[1])
    yres_mean=((Ymax-Ymin)/size[0])

    res=np.mean([xres_mean,yres_mean])
    
    lon_reg,lat_reg=np.meshgrid(np.arange(Xmin, Xmax, res),np.arange(Ymin, Ymax, res))
    
    logging.info('[InterpolatorCore] Source grid size : '+str(np.shape(sourceAxisX)))
    logging.info('[InterpolatorCore] Target grid size : '+str(np.shape(lon_reg)))

    return griddata((sourceAxisX.ravel(),sourceAxisY.ravel()), data.ravel(),(lon_reg, lat_reg), method=InterpolatorCore.HORIZONTAL_INTERPOLATION_METHOD)

def vertical_interpolation(sourceAxis,targetAxis,data):
    f = interp1d(sourceAxis,data,kind=InterpolatorCore.VERTICAL_INTERPOLATION_METHOD,bounds_error=False)
    return f(targetAxis)
