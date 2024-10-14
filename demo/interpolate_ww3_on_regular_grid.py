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
sys.path = ['../'] + sys.path

from spatialetl.coverage.TimeCoverage import TimeCoverage
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
     
    
    
    
       
        
    
    
    
    
    
