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
import time

from spatialetl.coverage.time_level_coverage import TimeLevelCoverage
from spatialetl.providers.symphonie.coverage.netcdf.v293.symphonie_reader import SYMPHONIEReader as CoverageReader
from spatialetl.providers.common.netcdf.coverage.default_writer import DefaultWriter
from spatialetl.utils.logger import logging

if __name__ == "__main__":
    start_time = time.time()
    logging.setLevel(logging.INFO)

    # Read file
    reader = CoverageReader('/work/sciences/projects/2022-Snale/simu-Aubenas/output/v2-1992/grid.nc',
                             '/work/sciences/projects/2022-Snale/simu-Aubenas/output/v2-1992/GRAPHIQUES/')

    #depths = [0.0, 10.0]
    coverage = TimeLevelCoverage(reader, resolution_x=0.001, resolution_y=0.001,resolution_z=0.0001,nb_thread=8);

    writer = DefaultWriter(coverage, '/tmp/symphonie_regular.nc')

    #writer.write_variable_baroclinic_sea_water_velocity()
    #writer.write_variable_barotropic_sea_water_velocity()
    #writer.write_variable_sea_water_column_thickness()
    #writer.write_variable_sea_surface_height_above_mean_sea_level()
    #writer.write_variable_wind_stress()
    #writer.write_variable_sea_water_temperature()
    #writer.write_variable_sea_water_salinity()
    # writer.write_variable_wet_binary_mask()
    writer.write_variable_2D_sea_binary_mask()
    #writer.write_variable_wind_10m()
    writer.write_variable_mesh_size()
    #writer.write_variable_sea_surface_wave_significant_height()
    #writer.write_variable_sea_surface_wave_mean_period()
    #writer.write_variable_sea_surface_wave_peak_period()
    #writer.write_variable_sea_surface_wave_from_direction()
    #writer.write_variable_sea_surface_wave_to_direction()
    writer.close()

    stop_time = time.time()
    print("----- Time ", (time.time() - start_time), " seconds -----")
    print('End of program')
     
    
    
    
       
        
    
    
    
    
    
