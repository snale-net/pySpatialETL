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
from datetime import datetime
import numpy as np


def time_interpolation(sourceAxis,targetAxis,data,method):
    logging.info("[InterpolatorCore] Interpolation for time : "+str(datetime.fromtimestamp(targetAxis[0]))+" with method '"+str(method)+"'.")
    f = interp1d(sourceAxis, data, kind=method, bounds_error=True)
    return f(targetAxis)
