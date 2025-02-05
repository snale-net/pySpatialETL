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

import numpy as np

#  Return the correlation coefficient of two time series (observed, modelized)
def correlation(obs,model):
  if len(obs) != len(model):
      raise ValueError("CORRELATION : obs and model timeseries have not the same length.")

  quot = np.sum( (obs-obs.mean()) * (model-model.mean()))
  div = np.sqrt ( np.sum((obs-obs.mean())**2) * np.sum((model-model.mean())**2) )
  return(quot/div)


# Return the relative bias of two time series (observed, simulated)
def bias(obs,model):
    if len(obs) != len(model):
      raise ValueError("BIAS : obs and model timeseries have not the same length.")

    return (model.mean()-obs.mean())/obs.mean()

# Return the root mean squared error of two time series (observed, simulated)
def rmse(obs,model):
    if len(obs) != len(model):
      raise ValueError("RMSE : obs and model timeseries have not the same length.")

    return (np.sqrt(np.sum((obs-model)**2)/len(obs)))

# Return the root mean squared relative error -- scatter index -- of two time series (observed, simulated)
# It presents the percentage of RMS difference with respect to mean observation or it gives the percentage of expected error for the parameter
def si(obs,model):
    if len(obs) != len(model):
      raise ValueError("SI : obs and model timeseries have not the same length.")

    return (np.sqrt( np.sum((obs-model)**2) / np.sum(obs*obs) ))

# Return the maximum error between two time series (observed, simulated)
def maxerr(obs,model):
    if len(obs) != len(model):
      raise ValueError("ERRMAX : obs and model timeseries have not the same length.")

    return (np.max(np.abs(obs-model)))
