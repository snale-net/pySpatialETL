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

# Lien vers le dossier de la lib
import sys
sys.path = ['/work/sciences/toolbox/python/pyGeoSpatialETL'] + sys.path

from coverage.TimeLevelCoverage import TimeLevelCoverage
from coverage.operator.interpolator.CoverageInterpolator import CoverageInterpolator
from coverage.operator.interpolator.InterpolatorCore import InterpolatorCore
from coverage.io.netcdf.symphonie.SymphonieReader2015 import SymphonieReader2015
from coverage.io.netcdf.DefaultWriter import DefaultWriter
import numpy as np
import logging

if __name__ == "__main__":
    print("Transform/Interpole Symphonie to GMT")
    
    logging.basicConfig(format='[%(levelname)s] %(message)s',level=logging.INFO)

    # Read file
    reader = SymphonieReader2015('/work/sciences/projects/WWB-2017/Manicouagan/configuration_V2015/TStra_N/OFFLINE/grid.nc',
                             '/work/sciences/projects/WWB-2017/Manicouagan/configuration_V2015/TStra_N/GRAPHIQUES/20090301_071217.nc')

    coverageOrig = TimeLevelCoverage(reader);

    # InterpolatorCore.INTERPOLATION_METHOD = "linear";
    InterpolatorCore.INTERPOLATION_METHOD = "nearest";

    depths = np.array([0.0, 10.0, 50.0, 100.0])
    coverage = TimeLevelCoverage(CoverageInterpolator(coverageOrig, 0.001, 0.001, depths))

    #Renommer toutes les noms des fonctions
    writer = DefaultWriter(coverage, '/tmp/symphonie_regular.nc')

    writer.write_variable_baroclinic_sea_water_velocity()
    writer.write_variable_barotropic_sea_water_velocity()
    writer.write_variable_sea_surface_height_above_mean_sea_level()
    writer.write_variable_wind_stress()
    writer.write_variable_sea_water_temperature()
    writer.write_variable_sea_water_salinity()
    # writer.write_variable_wet_binary_mask()
    writer.write_variable_2D_sea_binary_mask()
    writer.write_variable_wind_10m()
    writer.write_variable_mesh_size()
    #writer.write_variable_sea_surface_wave_significant_height()
    #writer.write_variable_sea_surface_wave_mean_period()
    #writer.write_variable_sea_surface_wave_peak_period()
    #writer.write_variable_sea_surface_wave_from_direction()
    #writer.write_variable_sea_surface_wave_to_direction()
    writer.close()
    
    print('End of programm')
     
    
    
    
       
        
    
    
    
    
    
