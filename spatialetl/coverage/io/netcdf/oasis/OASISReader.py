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
# Author : Fabien Rétif - fabien.retif@zoho.com
#
from __future__ import division, print_function, absolute_import
from spatialetl.coverage.io.CoverageReader import CoverageReader
from netCDF4 import Dataset, num2date
import numpy as np

class OASISReader(CoverageReader):
    """
La classe SymphonieReader permet de lire les données du format Symphonie

@param  myGrid : lien vers le fichier de grille (que l'on trouve dans le RDIR/tmp/grid.nc)
@param myFile : lien vers le fichier de données (d'après le namecouple)
"""
    def __init__(self,myGrid, myFile):
        CoverageReader.__init__(self,myFile);
        self.ncfile = Dataset(self.filename, 'r')
        self.grid = Dataset(myGrid, 'r')

    # Axis
    def read_axis_t(self,timestamp=0):
        data = self.ncfile.variables['time'][:]
        result = num2date(data, units = "seconds since 2013-03-05-18:00:00", calendar = "gregorian")

        if timestamp ==1:
            return [ (t - TimeCoverage.TIME_DATUM).total_seconds() \
                for t in result];
        else:
            return result

    def read_axis_x(self):
        return self.grid.variables['ww3t.lon'][:]

    def read_axis_y(self):
        return self.grid.variables['ww3t.lat'][:]

    # Data
    def read_variable_2D_mask(self):
        return self.grid.variables["mask_t"][0][:]

    def read_variable_bathymetry(self):
        return self.grid.variables["hm_w"][:]

    def read_variable_Ha(self):
        return self.ncfile.variables["Ha"][:]

    def read_variable_ssh_at_time(self,t):
        return self.ncfile.variables["ssh_w"][t][:]

    def read_variable_bathy_ssh_at_time(self,t):
        return self.ncfile.variables["hssh"][t][:]

    def read_variable_hs_at_time(self,t):
        return self.ncfile.variables["WW3__OHS"][t][:]

    def read_variable_waves_mean_period_at_time(self,t):
        return self.ncfile.variables["WW3_T0M1"][t][:]

    def read_variable_wetmask_at_time(self,t):
        return self.ncfile.variables["WW3_ODRY"][t][:]

    def read_variable_salinity_at_time_and_depth(self,index_t,index_z,depth,method="nearest"):
        mask_t = self.read_variable_3D_mask();
        lon_t = self.read_axis_x();
        lat_t = self.read_axis_y();
        xmax=np.shape(lon_t)[1]
        ymax=np.shape(lon_t)[0]
        data = self.ncfile.variables["sal"][index_t][:]
        sal = np.zeros([ymax,xmax])
        sal[:] = np.NAN

        for y in range(1,ymax-1):
            for x in range(1,xmax-1):

               if index_z[y,x] != -999 : # Le point (x,y) a une couche de profondeur depth

                    if mask_t[index_z[y,x],y,x] == 1.:
                        sal[y,x] = data[index_z[y,x],y,x]

        return sal


    def read_variable_surface_stokes_drift_at_time(self,index_t):

        lon_t = self.read_axis_x();
        lat_t = self.read_axis_y();
        data_u = self.ncfile.variables["WW3_USSX"][index_t][::]
        data_v = self.ncfile.variables["WW3_USSY"][index_t][::]

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

        # 1. On calcule les point à l'intérieur du domaine en excluant les bords
        for y in range(1,ymax-1):
            for x in range(1,xmax-1):

                # 1.1 On calcule la matrice de rotation
                x1=(lon_t[y,x+1]-lon_t[y,x-1])*np.pi/180.
                if(x1<-np.pi): x1=x1+2.*np.pi
                if(x1> np.pi): x1=x1-2.*np.pi
                x0=-np.arctan2((lat_t[y,x+1]-lat_t[y,x-1])*np.pi/180.,x1*np.cos(lat_t[y,x]*np.pi/180.))
                gridrotcos_t[y,x]=np.cos(x0)
                gridrotsin_t[y,x]=np.sin(x0)

                # u_left
                u_left = 0
                # Pas d'interpolation, on prend le plus proche inférieur
                u_left = data_u[y,x-1];

                # u_right
                u_right = 0
                # Pas d'interpolation, on prend le plus proche inférieur
                u_right = data_u[y,x];

                # v_down
                v_down = 0
                # Pas d'interpolation, on prend le plus proche inférieur
                v_down = data_v[y-1,x];

                # v_up
                v_up = 0
                # Pas d'interpolation, on prend le plus proche inférieur
                v_up = data_v[y,x];

                # 1.3 On calcule la demi-somme
                u_rot[y,x]=0.5*(u_left+u_right)
                v_rot[y,x]=0.5*(v_down+v_up)

                # 1.4 On applique la rotation
                #u_rot[y,x]=u[y,x]*gridrotcos_t[y,x]+v[y,x]*gridrotsin_t[y,x]
                #v_rot[y,x]=-u[y,x]*gridrotsin_t[y,x]+v[y,x]*gridrotcos_t[y,x]

        # 2. On duplique les point sur les bords.
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

    def read_variable_taw_at_time(self,index_t):

        lon_t = self.read_axis_x();
        lat_t = self.read_axis_y();
        data_u = self.ncfile.variables["WW3_TAWX"][index_t][::]
        data_v = self.ncfile.variables["WW3_TAWY"][index_t][::]

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

        # 1. On calcule les point à l'intérieur du domaine en excluant les bords
        for y in range(1,ymax-1):
            for x in range(1,xmax-1):

                # 1.1 On calcule la matrice de rotation
                x1=(lon_t[y,x+1]-lon_t[y,x-1])*np.pi/180.
                if(x1<-np.pi): x1=x1+2.*np.pi
                if(x1> np.pi): x1=x1-2.*np.pi
                x0=-np.arctan2((lat_t[y,x+1]-lat_t[y,x-1])*np.pi/180.,x1*np.cos(lat_t[y,x]*np.pi/180.))
                gridrotcos_t[y,x]=np.cos(x0)
                gridrotsin_t[y,x]=np.sin(x0)

                # u_left
                u_left = 0
                # Pas d'interpolation, on prend le plus proche inférieur
                u_left = data_u[y,x-1];

                # u_right
                u_right = 0
                # Pas d'interpolation, on prend le plus proche inférieur
                u_right = data_u[y,x];

                # v_down
                v_down = 0
                # Pas d'interpolation, on prend le plus proche inférieur
                v_down = data_v[y-1,x];

                # v_up
                v_up = 0
                # Pas d'interpolation, on prend le plus proche inférieur
                v_up = data_v[y,x];

                # 1.3 On calcule la demi-somme
                u_rot[y,x]=0.5*(u_left+u_right)
                v_rot[y,x]=0.5*(v_down+v_up)

                # 1.4 On applique la rotation
                #u_rot[y,x]=u[y,x]*gridrotcos_t[y,x]+v[y,x]*gridrotsin_t[y,x]
                #v_rot[y,x]=-u[y,x]*gridrotsin_t[y,x]+v[y,x]*gridrotcos_t[y,x]

        # 2. On duplique les point sur les bords.
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

    def read_variable_two_at_time(self,index_t):

        lon_t = self.read_axis_x();
        lat_t = self.read_axis_y();
        data_u = self.ncfile.variables["WW3_TWOX"][index_t][::]
        data_v = self.ncfile.variables["WW3_TWOY"][index_t][::]

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

        # 1. On calcule les point à l'intérieur du domaine en excluant les bords
        for y in range(1,ymax-1):
            for x in range(1,xmax-1):

                # 1.1 On calcule la matrice de rotation
                x1=(lon_t[y,x+1]-lon_t[y,x-1])*np.pi/180.
                if(x1<-np.pi): x1=x1+2.*np.pi
                if(x1> np.pi): x1=x1-2.*np.pi
                x0=-np.arctan2((lat_t[y,x+1]-lat_t[y,x-1])*np.pi/180.,x1*np.cos(lat_t[y,x]*np.pi/180.))
                gridrotcos_t[y,x]=np.cos(x0)
                gridrotsin_t[y,x]=np.sin(x0)

                # u_left
                u_left = 0
                # Pas d'interpolation, on prend le plus proche inférieur
                u_left = data_u[y,x-1];

                # u_right
                u_right = 0
                # Pas d'interpolation, on prend le plus proche inférieur
                u_right = data_u[y,x];

                # v_down
                v_down = 0
                # Pas d'interpolation, on prend le plus proche inférieur
                v_down = data_v[y-1,x];

                # v_up
                v_up = 0
                # Pas d'interpolation, on prend le plus proche inférieur
                v_up = data_v[y,x];

                # 1.3 On calcule la demi-somme
                u_rot[y,x]=0.5*(u_left+u_right)
                v_rot[y,x]=0.5*(v_down+v_up)

                # 1.4 On applique la rotation
                #u_rot[y,x]=u[y,x]*gridrotcos_t[y,x]+v[y,x]*gridrotsin_t[y,x]
                #v_rot[y,x]=-u[y,x]*gridrotsin_t[y,x]+v[y,x]*gridrotcos_t[y,x]

        # 2. On duplique les point sur les bords.
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

    def read_variable_waves_dir_at_time(self,index_t):

        lon_t = self.read_axis_x();
        lat_t = self.read_axis_y();
        data_u = self.ncfile.variables["WW3_SDIR"][index_t][::]
        data_v = self.ncfile.variables["WW3_CDIR"][index_t][::]

        xmax=np.shape(lon_t)[1]
        ymax=np.shape(lon_t)[0]
        gridrotcos_t = np.zeros([ymax,xmax])
        gridrotsin_t = np.zeros([ymax,xmax])

        u = np.zeros([ymax,xmax])
        u[:] = np.NAN
        v = np.zeros([ymax,xmax])
        v[:] = np.NAN
        rot = np.zeros([ymax,xmax])
        rot[:] = np.NAN

        # 1. On calcule les point à l'intérieur du domaine en excluant les bords
        for y in range(1,ymax-1):
            for x in range(1,xmax-1):

                # 1.1 On calcule la matrice de rotation
                x1=(lon_t[y,x+1]-lon_t[y,x-1])*np.pi/180.
                if(x1<-np.pi): x1=x1+2.*np.pi
                if(x1> np.pi): x1=x1-2.*np.pi
                x0=-np.arctan2((lat_t[y,x+1]-lat_t[y,x-1])*np.pi/180.,x1*np.cos(lat_t[y,x]*np.pi/180.))
                gridrotcos_t[y,x]=np.cos(x0)
                gridrotsin_t[y,x]=np.sin(x0)

                rad2deg = 180. / np.pi

                # 1.4 On applique la rotation
                #rot[y,x] = np.mod(270.- np.arctan2(data_u[y,x],data_v[y,x])+np.arctan2(gridrotsin_t[y,x],gridrotcos_t[y,x]),360.0)
                #rot[y,x] = np.mod(np.arctan2(data_u[y,x],data_v[y,x])+np.arctan2(gridrotsin_t[y,x],gridrotcos_t[y,x]),360.0)
                #rot[y,x] = 270. - np.arctan2(data_u[y,x],data_v[y,x])
                #rot[y,x] = np.mod ( 630. - rad2deg*np.arctan2(data_u[y,x],data_v[y,x])-np.arctan2(gridrotsin_t[y,x],gridrotcos_t[y,x])*rad2deg , 360. )
                #rot[y,x] = np.mod ( 630. - rad2deg*np.arctan2(data_u[y,x],data_v[y,x]), 360. )
                rot[y,x] = rad2deg*np.arctan2(data_u[y,x],data_v[y,x])

        # 2. On duplique les point sur les bords.
        # bottom
        rot[0,0:xmax]=rot[1,0:xmax]
        # up
        rot[ymax-1,0:xmax]=rot[ymax-2,0:xmax]
        # left
        rot[0:ymax,0]=rot[0:ymax,1]
        # right
        rot[0:ymax,xmax-1]=rot[0:ymax,xmax-2]

        return rot