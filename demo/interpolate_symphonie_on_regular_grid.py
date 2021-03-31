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

from spatialetl.coverage.TimeLevelCoverage import TimeLevelCoverage
from spatialetl.coverage.io.netcdf.symphonie.v293.SYMPHONIEReader import SYMPHONIEReader as CoverageReader
from spatialetl.coverage.io.netcdf.DefaultWriter import DefaultWriter
from spatialetl.utils.logger import logging

if __name__ == "__main__":
    logging.setLevel(logging.INFO)

    # Read file
    reader = CoverageReader('resources/symphonie_grid.nc',
                             'resources/symphonie_graphique.nc')

    coverageOrig = TimeLevelCoverage(reader);

    #depths = [0.0, 10.0]
    coverage = TimeLevelCoverage(reader, resolution_x=0.001, resolution_y=0.001,freq="3h",resolution_z=0.0001);

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
    
    print('End of program')
     
    
    
    
       
        
    
    
    
    
    
