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
# Author : Fabien Rétif - fabien.retif@zoho.com
#
from __future__ import division, print_function, absolute_import

import os

import numpy as np
from netCDF4 import Dataset, MFDataset, num2date

from spatialetl.coverage.TimeCoverage import TimeCoverage
from spatialetl.coverage.io.CoverageReader import CoverageReader


class MERCATORReader(CoverageReader):
    
    def __init__(self,m,d,t,u,v):
        CoverageReader.__init__(self,m)
        self.mask = Dataset(m, 'r')
        self.grid2D = Dataset(d, 'r')

        if os.path.isfile(t):
            self.gridT = Dataset(t, 'r')
        elif os.path.isdir(t):
            self.gridT = MFDataset(os.path.join(t,"*.nc"), 'r')
        elif t.endswith("*"):
            self.gridT = MFDataset(t+".nc", 'r')
        else:
            raise ValueError("Unable to decode file "+str(t))

        if os.path.isfile(u):
            self.gridU = Dataset(u, 'r')
        elif os.path.isdir(u):
            self.gridU = MFDataset(os.path.join(u,"*.nc"), 'r')
        elif u.endswith("*"):
            self.gridU = MFDataset(u+".nc", 'r')
        else:
            raise ValueError("Unable to decode file "+str(u))

        if os.path.isfile(v):
            self.gridV = Dataset(v, 'r')
        elif os.path.isdir(v):
            self.gridV = MFDataset(os.path.join(v,"*.nc"), 'r')
        elif v.endswith("*"):
            self.gridV = MFDataset(v+".nc", 'r')
        else:
            raise ValueError("Unable to decode file "+str(u))

    def close(self):
        self.mask.close()
        self.grid2D.close()
        self.gridT.close()
        self.gridU.close()
        self.gridV.close()

    def is_regular_grid(self):
        return False

    def get_x_size(self):
        return np.shape(self.gridT.variables['nav_lon'][:])[1];

    def get_y_size(self):
        return np.shape(self.gridT.variables['nav_lat'][:])[0];

    def get_z_size(self):
        return np.shape(self.gridT.variables['depth_t'][:])[0];

    def get_t_size(self):
        return np.shape(self.gridT.variables['time_counter'])[0];
     
    # Axis
    def read_axis_t(self,tmin,tmax,timestamp):
        """Attention si gridT, U,V,2D ont un time_counter different"""
        data = self.gridT.variables['time_counter'][tmin:tmax]
        result = num2date(data, units = self.gridT.variables['time_counter'].units, calendar = self.gridT.variables['time_counter'].calendar)
        
        if timestamp ==1:           
            return [ (t - TimeCoverage.TIME_DATUM).total_seconds() \
                for t in result];
        else:            
            return result
    
    def read_axis_x(self,xmin,xmax,ymin,ymax):
        return self.gridT.variables['nav_lon'][ymin:ymax,xmin:xmax]
    
    def read_axis_y(self,xmin,xmax,ymin,ymax):
        return self.gridT.variables['nav_lat'][ymin:ymax,xmin:xmax]
    
    def read_axis_z(self):       
        return self.gridT.variables['deptht']
    
    # Data    
    def read_variable_2D_sea_binary_mask(self,xmin,xmax,ymin,ymax):
        return self.mask.variables["tmask"][0][0][ymin:ymax,xmin:xmax]

    def read_variable_3D_sea_binary_mask(self,xmin,xmax,ymin,ymax):
        return self.mask.variables["tmask"][0][ymin:ymax,xmin:xmax]

    def read_variable_4D_sea_binary_mask(self,xmin,xmax,ymin,ymax):
        return self.mask.variables["tmask"][ymin:ymax,xmin:xmax]
    
    def read_variable_sea_surface_height_above_sea_level_at_time(self,t,xmin,xmax,ymin,ymax):
        return self.grid2D.variables["sossheig"][t,ymin:ymax,xmin:xmax]
     
    def read_variable_baroclinic_sea_water_velocity_at_time_and_depth(self,index_t,index_z,xmin,xmax,ymin,ymax):

        mask_t = self.read_variable_4D_sea_binary_mask();
        mask_u = self.mask.variables["umask"][:];
        mask_v = self.mask.variables["vmask"][:];
        lon_t = self.read_axis_x(xmin,xmax,ymin,ymax);
        lat_t = self.read_axis_y(xmin,xmax,ymin,ymax);
        data_u = self.gridU.variables["vozocrtx"][index_t,index_z,ymin:ymax,xmin:xmax]
        data_v = self.gridV.variables["vomecrty"][index_t,index_z,ymin:ymax,xmin:xmax]
        
        # compute and apply rotation matrix
        xmax=np.shape(lon_t)[1]
        ymax=np.shape(lon_t)[0]
        gridrotcos_t = np.zeros([ymax,xmax])
        gridrotsin_t = np.zeros([ymax,xmax])       
        
        u = np.zeros([ymax,xmax])
        u[:] = np.NAN
        v = np.zeros([ymax,xmax])
        v[:] = np.NAN
        u_rot = np.zeros([ymax,xmax])
        u_rot[:] = np.NAN
        v_rot = np.zeros([ymax,xmax])
        v_rot[:] = np.NAN

        # We process point inside the domain
        for y in range(1,ymax-1):
            for x in range(1,xmax-1):
                
                x1=(lon_t[y,x+1]-lon_t[y,x-1])*np.pi/180.
                if(x1<-np.pi): x1=x1+2.*np.pi
                if(x1> np.pi): x1=x1-2.*np.pi
                x0=-np.arctan2((lat_t[y,x+1]-lat_t[y,x-1])*np.pi/180.,x1*np.cos(lat_t[y,x]*np.pi/180.))
                gridrotcos_t[y,x]=np.cos(x0)
                gridrotsin_t[y,x]=np.sin(x0)

                if (mask_t[0,index_z[y,x],y,x] == 1.):
                    
                    u_left = 0
                    u_right = 0
                    v_down = 0
                    v_up = 0
                   
                    if (mask_u[0,index_z[y,x-1],y,x-1] == 1.):
                        u_left = data_u[index_z[y,x-1],y,x-1];

                    if (mask_u[0,index_z[y,x],y,x] == 1.):
                        u_right = data_u[index_z[y,x],y,x];

                    if (mask_v[0,index_z[y-1,x],y-1,x] == 1.):
                        v_down = data_v[index_z[y-1,x],y-1,x];

                    if (mask_v[0,index_z[y,x],y,x] == 1.):
                        v_up = data_v[index_z[y,x],y,x];

                    # compute an half-value
                    u[y,x]=0.5*(u_left+u_right)
                    v[y,x]=0.5*(v_down+v_up)
                    
                    # apply rotation                
                    u_rot[y,x]=u[y,x]*gridrotcos_t[y,x]+v[y,x]*gridrotsin_t[y,x]
                    v_rot[y,x]=-u[y,x]*gridrotsin_t[y,x]+v[y,x]*gridrotcos_t[y,x]   
          
        # We process boundaries point
        # bottom        
        u_rot[0,0:xmax]=u_rot[1,0:xmax]   
        v_rot[0,0:xmax]=v_rot[1,0:xmax] 
        # up
        u_rot[ymax-1,0:xmax]=u_rot[ymax-2,0:xmax] 
        v_rot[ymax-1,0:xmax]=v_rot[ymax-2,0:xmax] 
        
        # left
        u_rot[0:ymax,0]=u_rot[0:ymax,1]   
        v_rot[0:ymax,0]=v_rot[0:ymax,1]   
        # right
        u_rot[0:ymax,xmax-1]=u_rot[0:ymax,xmax-2]  
        v_rot[0:ymax,xmax-1]=v_rot[0:ymax,xmax-2]  

        return [u_rot,v_rot]

