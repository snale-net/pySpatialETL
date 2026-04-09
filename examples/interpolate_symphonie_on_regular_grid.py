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

import numpy

from spatialetl.coverage import Coverage
from spatialetl.coverage.time_level_coverage import TimeLevelCoverage
from spatialetl.providers.symphonie.coverage.netcdf.symphonie_reader import SYMPHONIEReader as CoverageReader
from spatialetl.providers.common.netcdf.coverage.default_writer import DefaultWriter
from spatialetl.utils.logger import logging

if __name__ == "__main__":
    start_time = time.time()
    logging.setLevel(logging.RUN)

    Coverage.HORIZONTAL_INTERPOLATOR = "cgal"

    # Read file
    #reader = CoverageReader('/work/sciences/projects/2022-Snale/simu-Aubenas/output/v2-1992/grid.nc',
    #                         '/work/sciences/projects/2022-Snale/simu-Aubenas/output/v2-1992/GRAPHIQUES/')
    reader = CoverageReader('/work/sciences/projects/Nokoue-2019/grid/grid_zoom.nc',
                            '/work/sciences/projects/Nokoue-2019/modelling/symphonie/outputs/REF2018/')

    depths = [0.0, 1.5]
    coverage = TimeLevelCoverage(reader, resolution_x=0.0005, resolution_y=0.0005,resolution_z=0.5, zbox=depths,freq="1d",nb_thread=15);

    writer = DefaultWriter(coverage, '/work/sciences/projects/Nokoue-2019/modelling/symphonie/outputs/symphonie_regular.nc')

    #writer.write_variable_baroclinic_sea_water_velocity()
    writer.write_variable_barotropic_sea_water_velocity()
    #writer.write_variable_sea_water_column_thickness()
    writer.write_variable_sea_surface_height_above_mean_sea_level()
    #writer.write_variable_wind_stress()
    #writer.write_variable_sea_water_temperature()
    #writer.write_variable_sea_water_salinity()
    writer.write_variable_sea_surface_salinity()
    writer.write_variable_sea_water_salinity_at_ground_level()
    # writer.write_variable_wet_binary_mask()
    #writer.write_variable_2D_sea_binary_mask()
    #writer.write_variable_wind_10m()
    #writer.write_variable_mesh_size()
    #writer.write_variable_bathymetry()
    #writer.write_variable_sea_surface_wave_significant_height()
    #writer.write_variable_sea_surface_wave_mean_period()
    #writer.write_variable_sea_surface_wave_peak_period()
    #writer.write_variable_sea_surface_wave_from_direction()
    #writer.write_variable_sea_surface_wave_to_direction()
    writer.close()

    stop_time = time.time()
    print("----- Time ", (time.time() - start_time), " seconds -----")
    print('End of program')
     
    
    
    
       
        
    
    
    
    
    
