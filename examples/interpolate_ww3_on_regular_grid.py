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
from spatialetl.coverage.time_coverage import TimeCoverage
from spatialetl.coverage.operator.interpolator.CoverageInterpolator import CoverageInterpolator
from spatialetl.coverage.operator.interpolator.InterpolatorCore import InterpolatorCore
from spatialetl.coverage.io.netcdf.ww3.WW3Reader import WW3Reader
from spatialetl.coverage.io.netcdf.DefaultWriter import DefaultWriter
import numpy as np
import logging

if __name__ == "__main__":
    print("Transform/Interpole Symphonie to GMT")
    
    logging.basicConfig(format='[%(levelname)s] %(message)s',level=logging.INFO)

    # Read file
    reader = WW3Reader('/work/sciences/projects/ArcEmeraude_2018/waves/Re_ Bon___/ww3.20150105_T405_2s.nc')

    coverageOrig = TimeCoverage(reader);

    # InterpolatorCore.INTERPOLATION_METHOD = "linear";
    InterpolatorCore.INTERPOLATION_METHOD = "nearest";

    depths = np.array([])
    coverage = TimeCoverage(CoverageInterpolator(coverageOrig, 0.1, 0.1, depths))

    #Renommer toutes les noms des fonctions
    writer = DefaultWriter(coverage, '/tmp/ww3_regular.nc')

    writer.write_variable_barotropic_sea_water_velocity()
    writer.write_variable_sea_surface_height_above_mean_sea_level()
    writer.write_variable_2D_sea_binary_mask()
    writer.write_variable_sea_surface_wave_significant_height()
    writer.write_variable_sea_surface_wave_mean_period()
    # writer.write_variable_sea_surface_wave_peak_period()
    writer.write_variable_sea_surface_wave_from_direction()
    writer.write_variable_sea_surface_wave_to_direction()
    writer.write_variable_atmosphere_momentum_flux_to_waves()
    writer.write_variable_waves_momentum_flux_to_ocean()
    writer.write_variable_sea_surface_wave_stokes_drift_velocity()
    writer.write_variable_radiation_pressure_bernouilli_head()
    writer.write_variable_wind_10m()
    #writer.write_variable_sea_surface_wave_energy_dissipation_at_ground_level()
    writer.write_variable_sea_surface_wave_breaking_height()
    writer.write_variable_sea_surface_wave_energy_flux_to_ocean()
    writer.close()
    
    print('End of program')
     
    
    
    
       
        
    
    
    
    
    
