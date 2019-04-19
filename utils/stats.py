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
