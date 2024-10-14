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
from spatialetl.coverage.io.netcdf.mercator.v2020.MERCATORReader import MERCATORReader
from spatialetl.coverage.io.netcdf.ww3.WW3Writer import WW3Writer
import logging
import numpy as np

if __name__ == "__main__":
    """
    Cette routine permet de transformer les données du format Mercator à WaveWatch III en vu de forcer le modèle de vagues.
    """
    print("Transform MERCATOR to WW3")
    
    logging.basicConfig(format='[%(levelname)s] %(message)s',level=logging.INFO)
    
    # Read file
    reader = MERCATORReader('/home/fabien/mercator/ext-PSY2V4R4_mask.nc',
                            '/home/fabien/mercator/ext-PSY2V4R4_1dAV_20130201_20130202_grid2D_R20130213.nc',
                            '/home/fabien/mercator/ext-PSY2V4R4_1dAV_20130201_20130202_gridT_R20130213.nc',
                            '/home/fabien/mercator/ext-PSY2V4R4_1dAV_20130201_20130202_gridU_R20130213.nc',
                            '/home/fabien/mercator/ext-PSY2V4R4_1dAV_20130201_20130202_gridV_R20130213.nc')
        
    coverage = TimeLevelCoverage(reader);
    depths = np.array([0.0, 10.0, 50.0, 100.0])
    
    writer = WW3Writer(coverage,'/tmp/ww3.mercator-forcing.nc',depths)

    writer.write_variable_sea_surface_height_above_mean_sea_level()
    writer.write_variable_barotropic_sea_water_velocity();
    writer.close()    
    
    print('End of program')
     
    
    
    
       
        
    
    
    
    
    
