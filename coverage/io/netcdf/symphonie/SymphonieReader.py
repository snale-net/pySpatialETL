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
from coverage.io.CoverageReader import CoverageReader
from coverage.TimeCoverage import TimeCoverage
from netCDF4 import Dataset, MFDataset, num2date
import numpy as np
import os

class SymphonieReader(CoverageReader):
    """
La classe SymphonieReader permet de lire les données du format Symphonie

@param  myGrid : lien vers le fichier de grille (que l'on trouve dans le RDIR/tmp/grid.nc)
@param myFile : lien vers le fichier de données (que l'on trouve dans GRAPHIQUES)
"""
    def __init__(self,myGrid, myFile):   
        CoverageReader.__init__(self,myFile);

        if os.path.isfile(self.filename):
            self.ncfile = Dataset(self.filename, 'r')
        elif os.path.isdir(self.filename):
            self.ncfile = MFDataset(os.path.join(self.filename,"*.nc"), 'r')

        self.grid = Dataset(myGrid, 'r')
        lon_t = self.read_axis_x()
        lat_t = self.read_axis_y()
        self.xmax = np.shape(lon_t)[1]
        self.ymax = np.shape(lon_t)[0]
        self.zmax = np.shape(self.read_axis_z())[0]
        self.gridrotcos_t = np.zeros([self.ymax, self.xmax])
        self.gridrotsin_t = np.zeros([self.ymax, self.xmax])

        for y in range(1, self.ymax - 1):
            for x in range(1, self.xmax - 1):

                x1 = (lon_t[y, x + 1] - lon_t[y, x - 1]) * np.pi / 180.
                if (x1 < -np.pi): x1 = x1 + 2. * np.pi
                if (x1 > np.pi): x1 = x1 - 2. * np.pi
                x0 = -np.arctan2((lat_t[y, x + 1] - lat_t[y, x - 1]) * np.pi / 180.,
                                 x1 * np.cos(lat_t[y, x] * np.pi / 180.))
                self.gridrotcos_t[y, x] = np.cos(x0)
                self.gridrotsin_t[y, x] = np.sin(x0)

    def read_axis_x(self):
        return self.grid.variables['longitude_t'][:]

    def read_axis_y(self):
        return self.grid.variables['latitude_t'][:]

    def read_axis_z(self):
        lev = np.ma.masked_where(self.grid.variables["mask_t"][::] < 1 , self.grid.variables['depth_t'][::])
        np.ma.set_fill_value(lev,np.nan)
        lev[::] *= -1.0  # inverse la profondeur
        return lev.filled()

    def read_axis_t(self, timestamp):
        data = self.ncfile.variables['time'][:]
        result = num2date(data, units=self.ncfile.variables['time'].units.replace('from', 'since').replace('jan',
                                                                                                           '01').replace(
            'feb', '02').replace('mar', '03').replace('apr', '04').replace('may', '05').replace('jun', '06').replace(
            'jul', '07').replace('aug', '08').replace('sept', '09').replace('oct', '10').replace('nov', '11').replace(
            'dec', '12'), calendar=self.ncfile.variables['time'].calendar)

        if timestamp == 1:
            return [(t - TimeCoverage.TIME_DATUM).total_seconds() \
                    for t in result];
        else:
            return result

    # Variables
    def read_variable_time(self):
        return self.read_axis_t(timestamp=0)

    def read_variable_2D_sea_binary_mask(self):
        index_z = self.zmax - 1
        return self.grid.variables["mask_t"][index_z][:]
        #data = np.ma.masked_where(self.grid.variables["mask_t"][index_z][:] < 1,
        #                          self.grid.variables["mask_t"][index_z][:])
        #np.ma.set_fill_value(data, -9999)
        #return data.filled()

    def read_variable_3D_sea_binary_mask(self):
        return self.grid.variables["mask_t"][::]
        #data = np.ma.masked_where(self.grid.variables["mask_t"][::] < 1,
        #                          self.grid.variables["mask_t"][::])
        #np.ma.set_fill_value(data, -9999)
        #return data.filled()

    def read_variable_wet_binary_mask_at_time(self, index_t):
        return self.ncfile.variables["wetmask_t"][index_t][:]

    def read_variable_mesh_size(self):
        data = self.grid.variables["sqrt_dxdy"][:]
        data[data < 0] = np.nan
        return data

    #################
    # HYDRO
    # Sea Surface
    #################

    def read_variable_sea_surface_height_above_mean_sea_level_at_time(self, index_t):
        index_z = self.zmax - 1
        data = np.ma.masked_where(self.grid.variables["mask_t"][index_z][:] < 1 ,self.ncfile.variables["ssh_w"][index_t][:])
        np.ma.set_fill_value(data, np.nan)
        return data.filled()

    def read_variable_sea_surface_temperature_at_time(self, index_t):
        index_z = self.zmax - 1
        #return self.ncfile.variables["tem"][index_t][index_z][:]
        data = np.ma.masked_where(self.grid.variables["mask_t"][index_z][:] < 1,
                                  self.ncfile.variables["tem"][index_t][index_z][:])
        np.ma.set_fill_value(data, np.nan)
        return data.filled()


        # mask_t = self.read_variable_3D_sea_binary_mask();
        # lon_t = self.read_axis_x();
        # xmax = np.shape(lon_t)[1]
        # ymax = np.shape(lon_t)[0]
        # data = self.ncfile.variables["tem"][index_t][:]
        # result = np.zeros([ymax, xmax])
        # result[:] = np.NAN
        #
        # for y in range(0, ymax):
        #     for x in range(0, xmax):
        #
        #         if mask_t[index_z, y, x] == 1.:
        #             result[y, x] = data[index_z, y, x]
        #
        # return result

    def read_variable_sea_surface_salinity_at_time(self, index_t):
        index_z = self.zmax - 1
        data = np.ma.masked_where(self.grid.variables["mask_t"][index_z][:] < 1,
                                  self.ncfile.variables["sal"][index_t][index_z][:])
        np.ma.set_fill_value(data, np.nan)
        return data.filled()

        # mask_t = self.read_variable_3D_sea_binary_mask();
        # lon_t = self.read_axis_x();
        # xmax = np.shape(lon_t)[1]
        # ymax = np.shape(lon_t)[0]
        # data = self.ncfile.variables["sal"][index_t][:]
        # result = np.zeros([ymax, xmax])
        # result[:] = np.NAN
        #
        # for y in range(0, ymax):
        #     for x in range(0, xmax):
        #
        #         if mask_t[index_z, y, x] == 1.:
        #             result[y, x] = data[index_z, y, x]
        #
        # return result


    def read_variable_sea_water_velocity_at_sea_water_surface_at_time(self, index_t):
        index_z = self.zmax - 1
        mask_t = self.read_variable_3D_sea_binary_mask();
        mask_u = self.grid.variables["mask_u"][index_z][:];
        mask_v = self.grid.variables["mask_v"][index_z][:];
        data_u = self.ncfile.variables["vel_u"][index_t][index_z][:]
        data_v = self.ncfile.variables["vel_v"][index_t][index_z][:]

        u = np.zeros([self.ymax, self.xmax])
        u[:] = np.NAN
        v = np.zeros([self.ymax, self.xmax])
        v[:] = np.NAN
        u_rot = np.zeros([self.ymax, self.xmax])
        u_rot[:] = np.NAN
        v_rot = np.zeros([self.ymax, self.xmax])
        v_rot[:] = np.NAN

        # 1. On calcule les point à l'intérieur du domaine en excluant les bords
        for y in range(1, self.ymax - 1):
            for x in range(1, self.xmax - 1):

                if mask_t[index_z, y, x] == 1.:

                    # 1.2 On récupère les valeurs aux point encadrant X pour faire la demi-somme
                    ##############################
                    #           v_up
                    #
                    #   u_left   X     u_right
                    #
                    #         v_bottom
                    #############################

                    # u_left
                    u_left = 0
                    if mask_u[y, x - 1] == 1.:
                        u_left = data_u[y, x - 1];

                    # u_right
                    u_right = 0
                    if mask_u[y, x] == 1.:
                        u_right = data_u[y, x];

                    # v_down
                    v_down = 0
                    if mask_v[y - 1, x] == 1.:
                        v_down = data_v[y - 1, x];

                    # v_up
                    v_up = 0
                    if mask_v[y, x] == 1.:
                        v_up = data_v[y, x];

                    # 1.3 On calcule la demi-somme
                    u[y, x] = 0.5 * (u_left + u_right)
                    v[y, x] = 0.5 * (v_down + v_up)

                    # 1.4 On applique la rotation
                    u_rot[y, x] = u[y, x] * self.gridrotcos_t[y, x] + v[y, x] * self.gridrotsin_t[y, x]
                    v_rot[y, x] = -u[y, x] * self.gridrotsin_t[y, x] + v[y, x] * self.gridrotcos_t[y, x]

        # 2. On duplique les point sur les bords.
        # bottom
        u_rot[0, 0:self.xmax] = u_rot[1, 0:self.xmax]
        v_rot[0, 0:self.xmax] = v_rot[1, 0:self.xmax]
        # up
        u_rot[self.ymax - 1, 0:self.xmax] = u_rot[self.ymax - 2, 0:self.xmax]
        v_rot[self.ymax - 1, 0:self.xmax] = v_rot[self.ymax - 2, 0:self.xmax]

        # left
        u_rot[0:self.ymax, 0] = u_rot[0:self.ymax, 1]
        v_rot[0:self.ymax, 0] = v_rot[0:self.ymax, 1]
        # right
        u_rot[0:self.ymax, self.xmax - 1] = u_rot[0:self.ymax, self.xmax - 2]
        v_rot[0:self.ymax, self.xmax - 1] = v_rot[0:self.ymax, self.xmax - 2]

        return [u_rot, v_rot]


    #################
    # HYDRO
    # Ground level
    #################

    def read_variable_sea_water_temperature_at_ground_level_at_time(self, index_t):
        index_z = 0
        data = np.ma.masked_where(self.grid.variables["mask_t"][index_z][:] < 1,
                                  self.ncfile.variables["tem"][index_t][index_z][:])
        np.ma.set_fill_value(data, np.nan)
        return data.filled()

        # mask_t = self.read_variable_3D_sea_binary_mask();
        # lon_t = self.read_axis_x();
        # xmax = np.shape(lon_t)[1]
        # ymax = np.shape(lon_t)[0]
        # data = self.ncfile.variables["tem"][index_t][:]
        # result = np.zeros([ymax, xmax])
        # result[:] = np.NAN
        #
        #
        # for y in range(0, ymax):
        #     for x in range(0, xmax):
        #
        #         if mask_t[index_z, y, x] == 1.:
        #             result[y, x] = data[index_z, y, x]
        #
        # return result

    def read_variable_sea_water_salinity_at_ground_level_at_time(self, index_t):
        index_z = 0
        data = np.ma.masked_where(self.grid.variables["mask_t"][index_z][:] < 1,
                                  self.ncfile.variables["sal"][index_t][index_z][:])
        np.ma.set_fill_value(data, np.nan)
        return data.filled()

        # mask_t = self.read_variable_3D_sea_binary_mask();
        # lon_t = self.read_axis_x();
        # xmax = np.shape(lon_t)[1]
        # ymax = np.shape(lon_t)[0]
        # data = self.ncfile.variables["sal"][index_t][:]
        # result = np.zeros([ymax, xmax])
        # result[:] = np.NAN
        #
        #
        # for y in range(0, ymax):
        #     for x in range(0, xmax):
        #
        #         if mask_t[index_z, y, x] == 1.:
        #             result[y, x] = data[index_z, y, x]
        #
        # return result


    def read_variable_sea_water_velocity_at_ground_level_at_time(self, index_t):
        index_z = 0
        mask_t = self.read_variable_3D_sea_binary_mask();
        mask_u = self.grid.variables["mask_u"][index_z][:];
        mask_v = self.grid.variables["mask_v"][index_z][:];
        data_u = self.ncfile.variables["vel_u"][index_t][index_z][:]
        data_v = self.ncfile.variables["vel_v"][index_t][index_z][:]

        u = np.zeros([self.ymax, self.xmax])
        u[:] = np.NAN
        v = np.zeros([self.ymax, self.xmax])
        v[:] = np.NAN
        u_rot = np.zeros([self.ymax, self.xmax])
        u_rot[:] = np.NAN
        v_rot = np.zeros([self.ymax, self.xmax])
        v_rot[:] = np.NAN

        # 1. On calcule les point à l'intérieur du domaine en excluant les bords
        for y in range(1, self.ymax - 1):
            for x in range(1, self.xmax - 1):

                if mask_t[index_z, y, x] == 1.:

                    # 1.2 On récupère les valeurs aux point encadrant X pour faire la demi-somme
                    ##############################
                    #           v_up
                    #
                    #   u_left   X     u_right
                    #
                    #         v_bottom
                    #############################

                    # u_left
                    u_left = 0
                    if mask_u[y, x - 1] == 1.:
                        u_left = data_u[y, x - 1];

                    # u_right
                    u_right = 0
                    if mask_u[y, x] == 1.:
                        u_right = data_u[y, x];

                    # v_down
                    v_down = 0
                    if mask_v[y - 1, x] == 1.:
                        v_down = data_v[y - 1, x];

                    # v_up
                    v_up = 0
                    if mask_v[y, x] == 1.:
                        v_up = data_v[y, x];

                    # 1.3 On calcule la demi-somme
                    u[y, x] = 0.5 * (u_left + u_right)
                    v[y, x] = 0.5 * (v_down + v_up)

                    # 1.4 On applique la rotation
                    u_rot[y, x] = u[y, x] * self.gridrotcos_t[y, x] + v[y, x] * self.gridrotsin_t[y, x]
                    v_rot[y, x] = -u[y, x] * self.gridrotsin_t[y, x] + v[y, x] * self.gridrotcos_t[y, x]

        # 2. On duplique les point sur les bords.
        # bottom
        u_rot[0, 0:self.xmax] = u_rot[1, 0:self.xmax]
        v_rot[0, 0:self.xmax] = v_rot[1, 0:self.xmax]
        # up
        u_rot[self.ymax - 1, 0:self.xmax] = u_rot[self.ymax - 2, 0:self.xmax]
        v_rot[self.ymax - 1, 0:self.xmax] = v_rot[self.ymax - 2, 0:self.xmax]

        # left
        u_rot[0:self.ymax, 0] = u_rot[0:self.ymax, 1]
        v_rot[0:self.ymax, 0] = v_rot[0:self.ymax, 1]
        # right
        u_rot[0:self.ymax, self.xmax - 1] = u_rot[0:self.ymax, self.xmax - 2]
        v_rot[0:self.ymax, self.xmax - 1] = v_rot[0:self.ymax, self.xmax - 2]

        return [u_rot, v_rot]

    #################
    # HYDRO
    # 2D
    #################

    def read_variable_bathymetry(self):
        return self.grid.variables["hm_w"][:]

    def read_variable_barotropic_sea_water_velocity_at_time(self, index_t):
        index_z = self.zmax - 1
        mask_t = self.read_variable_2D_sea_binary_mask();
        mask_u = self.grid.variables["mask_u"][index_z][:];
        mask_v = self.grid.variables["mask_v"][index_z][:];
        data_u = self.ncfile.variables["velbar_u"][index_t][::]
        data_v = self.ncfile.variables["velbar_v"][index_t][::]
        u = np.zeros([self.ymax, self.xmax])
        u[:] = np.NAN
        v = np.zeros([self.ymax, self.xmax])
        v[:] = np.NAN
        u_rot = np.zeros([self.ymax, self.xmax])
        u_rot[:] = np.NAN
        v_rot = np.zeros([self.ymax, self.xmax])
        v_rot[:] = np.NAN

        # 1. On calcule les point à l'intérieur du domaine en excluant les bords
        for y in range(1, self.ymax - 1):
            for x in range(1, self.xmax - 1):

                if mask_t[y, x] == 1.:

                    # 1.2 On récupère les valeurs aux point encadrant X pour faire la demi-somme
                    ##############################
                    #           v_up
                    #
                    #   u_left   X     u_right
                    #
                    #         v_bottom
                    #############################

                    # u_left
                    u_left = 0
                    if mask_u[y, x - 1] == 1.:
                        u_left = data_u[y, x - 1];

                    # u_right
                    u_right = 0
                    if mask_u[y, x] == 1.:
                        u_right = data_u[y, x];

                    # v_down
                    v_down = 0
                    if mask_v[y - 1, x] == 1.:
                        v_down = data_v[y - 1, x];

                    # v_up
                    v_up = 0
                    if mask_v[y, x] == 1.:
                        v_up = data_v[y, x];

                    # 1.3 On calcule la demi-somme
                    u[y, x] = 0.5 * (u_left + u_right)
                    v[y, x] = 0.5 * (v_down + v_up)

                    # 1.4 On applique la rotation
                    u_rot[y, x] = u[y, x] * self.gridrotcos_t[y, x] + v[y, x] * self.gridrotsin_t[y, x]
                    v_rot[y, x] = -u[y, x] * self.gridrotsin_t[y, x] + v[y, x] * self.gridrotcos_t[y, x]

        # 2. On duplique les point sur les bords.
        # bottom
        u_rot[0, 0:self.xmax] = u_rot[1, 0:self.xmax]
        v_rot[0, 0:self.xmax] = v_rot[1, 0:self.xmax]
        # up
        u_rot[self.ymax - 1, 0:self.xmax] = u_rot[self.ymax - 2, 0:self.xmax]
        v_rot[self.ymax - 1, 0:self.xmax] = v_rot[self.ymax - 2, 0:self.xmax]

        # left
        u_rot[0:self.ymax, 0] = u_rot[0:self.ymax, 1]
        v_rot[0:self.ymax, 0] = v_rot[0:self.ymax, 1]
        # right
        u_rot[0:self.ymax, self.xmax - 1] = u_rot[0:self.ymax, self.xmax - 2]
        v_rot[0:self.ymax, self.xmax - 1] = v_rot[0:self.ymax, self.xmax - 2]

        return [u_rot, v_rot]

    #################
    # HYDRO
    # 3D
    #################
    def read_variable_sea_water_temperature_at_time_and_depth(self, index_t, index_z):
        data = np.ma.masked_where(self.grid.variables["mask_t"][index_z][:] < 1,
                                  self.ncfile.variables["tem"][index_t][index_z][:])
        np.ma.set_fill_value(data, np.nan)
        return data.filled()

    def read_variable_sea_water_salinity_at_time_and_depth(self, index_t, index_z):
        data = np.ma.masked_where(self.grid.variables["mask_t"][index_z][:] < 1,
                                  self.ncfile.variables["sal"][index_t][index_z][:])
        np.ma.set_fill_value(data, np.nan)
        return data.filled()

    def read_variable_baroclinic_sea_water_velocity_at_time_and_depth(self, index_t, index_z):
        mask_t = self.read_variable_3D_sea_binary_mask();
        mask_u = self.grid.variables["mask_u"][index_z][:];
        mask_v = self.grid.variables["mask_v"][index_z][:];
        data_u = self.ncfile.variables["vel_u"][index_t][index_z][:]
        data_v = self.ncfile.variables["vel_v"][index_t][index_z][:]

        u = np.zeros([self.ymax, self.xmax])
        u[:] = np.NAN
        v = np.zeros([self.ymax, self.xmax])
        v[:] = np.NAN
        u_rot = np.zeros([self.ymax, self.xmax])
        u_rot[:] = np.NAN
        v_rot = np.zeros([self.ymax, self.xmax])
        v_rot[:] = np.NAN

        # 1. On calcule les point à l'intérieur du domaine en excluant les bords
        for y in range(1, self.ymax - 1):
            for x in range(1, self.xmax - 1):

                if mask_t[index_z, y, x] == 1.:

                    # 1.2 On récupère les valeurs aux point encadrant X pour faire la demi-somme
                    ##############################
                    #           v_up
                    #
                    #   u_left   X     u_right
                    #
                    #         v_bottom
                    #############################

                    # u_left
                    u_left = 0
                    if mask_u[y, x - 1] == 1.:
                        u_left = data_u[y, x - 1];

                    # u_right
                    u_right = 0
                    if mask_u[y, x] == 1.:
                        u_right = data_u[y, x];

                    # v_down
                    v_down = 0
                    if mask_v[y - 1, x] == 1.:
                        v_down = data_v[y - 1, x];

                    # v_up
                    v_up = 0
                    if mask_v[y, x] == 1.:
                        v_up = data_v[y, x];

                    # 1.3 On calcule la demi-somme
                    u[y, x] = 0.5 * (u_left + u_right)
                    v[y, x] = 0.5 * (v_down + v_up)

                    # 1.4 On applique la rotation
                    u_rot[y, x] = u[y, x] * self.gridrotcos_t[y, x] + v[y, x] * self.gridrotsin_t[y, x]
                    v_rot[y, x] = -u[y, x] * self.gridrotsin_t[y, x] + v[y, x] * self.gridrotcos_t[y, x]

        # 2. On duplique les point sur les bords.
        # bottom
        u_rot[0, 0:self.xmax] = u_rot[1, 0:self.xmax]
        v_rot[0, 0:self.xmax] = v_rot[1, 0:self.xmax]
        # up
        u_rot[self.ymax - 1, 0:self.xmax] = u_rot[self.ymax - 2, 0:self.xmax]
        v_rot[self.ymax - 1, 0:self.xmax] = v_rot[self.ymax - 2, 0:self.xmax]

        # left
        u_rot[0:self.ymax, 0] = u_rot[0:self.ymax, 1]
        v_rot[0:self.ymax, 0] = v_rot[0:self.ymax, 1]
        # right
        u_rot[0:self.ymax, self.xmax - 1] = u_rot[0:self.ymax, self.xmax - 2]
        v_rot[0:self.ymax, self.xmax - 1] = v_rot[0:self.ymax, self.xmax - 2]

        return [u_rot, v_rot]

    #################
    # WAVES
    # Sea Surface
    #################

    def read_variable_sea_surface_wave_significant_height_at_time(self, index_t):
        return self.ncfile.variables["hs_wave_t"][index_t][:]

    def read_variable_sea_surface_wave_mean_period_at_time(self, index_t):
        return self.ncfile.variables["t_wave_t"][index_t][:]

    def read_variable_sea_surface_wave_to_direction_at_time(self, index_t):
        return self.ncfile.variables["dir_wave_t"][index_t][:]

    def read_variable_sea_surface_wave_stokes_drift_velocity_at_time(self, index_t):
        index_z = self.zmax - 1
        mask_t = self.read_variable_2D_sea_binary_mask();
        mask_u = self.grid.variables["mask_u"][index_z][:];
        mask_v = self.grid.variables["mask_v"][index_z][:];
        data_u = self.ncfile.variables["velbarstokes_u"][index_t][::]
        data_v = self.ncfile.variables["velbarstokes_v"][index_t][::]

        u = np.zeros([self.ymax, self.xmax])
        u[:] = np.NAN
        v = np.zeros([self.ymax, self.xmax])
        v[:] = np.NAN
        u_rot = np.zeros([self.ymax, self.xmax])
        u_rot[:] = np.NAN
        v_rot = np.zeros([self.ymax, self.xmax])
        v_rot[:] = np.NAN

        # 1. On calcule les point à l'intérieur du domaine en excluant les bords
        for y in range(1, self.ymax - 1):
            for x in range(1, self.xmax - 1):

                if mask_t[y, x] == 1.:

                    # 1.2 On récupère les valeurs aux point encadrant X pour faire la demi-somme.
                    ##############################
                    #           v_up
                    #
                    #   u_left   X     u_right
                    #
                    #         v_bottom
                    #############################

                    # u_left
                    u_left = 0
                    if mask_u[y, x - 1] == 1.:
                        u_left = data_u[y, x - 1];

                    # u_right
                    u_right = 0
                    if mask_u[y, x] == 1.:
                        u_right = data_u[y, x];

                    # v_down
                    v_down = 0
                    if mask_v[y - 1, x] == 1.:
                        v_down = data_v[y - 1, x];

                    # v_up
                    v_up = 0
                    if mask_v[y, x] == 1.:
                        v_up = data_v[y, x];

                    # 1.3 On calcule la demi-somme
                    u[y, x] = 0.5 * (u_left + u_right)
                    v[y, x] = 0.5 * (v_down + v_up)

                    # 1.4 On applique la rotation
                    u_rot[y, x] = u[y, x] * self.gridrotcos_t[y, x] + v[y, x] * self.gridrotsin_t[y, x]
                    v_rot[y, x] = -u[y, x] * self.gridrotsin_t[y, x] + v[y, x] * self.gridrotcos_t[y, x]

        # 2. On duplique les point sur les bords.
        # bottom
        u_rot[0, 0:self.xmax] = u_rot[1, 0:self.xmax]
        v_rot[0, 0:self.xmax] = v_rot[1, 0:self.xmax]
        # up
        u_rot[self.ymax - 1, 0:self.xmax] = u_rot[self.ymax - 2, 0:self.xmax]
        v_rot[self.ymax - 1, 0:self.xmax] = v_rot[self.ymax - 2, 0:self.xmax]

        # left
        u_rot[0:self.ymax, 0] = u_rot[0:self.ymax, 1]
        v_rot[0:self.ymax, 0] = v_rot[0:self.ymax, 1]
        # right
        u_rot[0:self.ymax, self.xmax - 1] = u_rot[0:self.ymax, self.xmax - 2]
        v_rot[0:self.ymax, self.xmax - 1] = v_rot[0:self.ymax, self.xmax - 2]

        return [u_rot, v_rot]


    #################
    # WAVES
    # Momentum flux
    #################

    def read_variable_atmosphere_momentum_flux_to_waves_at_time(self, index_t):
        index_z = self.zmax - 1
        mask_t = self.read_variable_2D_sea_binary_mask();
        mask_u = self.grid.variables["mask_u"][index_z][:];
        mask_v = self.grid.variables["mask_v"][index_z][:];
        data_u = self.ncfile.variables["tawx"][index_t][::]
        data_v = self.ncfile.variables["tawy"][index_t][::]

        u = np.zeros([self.ymax, self.xmax])
        u[:] = np.NAN
        v = np.zeros([self.ymax, self.xmax])
        v[:] = np.NAN
        u_rot = np.zeros([self.ymax, self.xmax])
        u_rot[:] = np.NAN
        v_rot = np.zeros([self.ymax, self.xmax])
        v_rot[:] = np.NAN

        # 1. On calcule les point à l'intérieur du domaine en excluant les bords
        for y in range(1, self.ymax - 1):
            for x in range(1, self.xmax - 1):

                if mask_t[y, x] == 1.:

                    # 1.2 On récupère les valeurs aux point encadrant X pour faire la demi-somme
                    ##############################
                    #           v_up
                    #
                    #   u_left   X     u_right
                    #
                    #         v_bottom
                    #############################

                    # u_left
                    u_left = 0
                    if mask_u[y, x - 1] == 1.:
                        u_left = data_u[y, x - 1];

                    # u_right
                    u_right = 0
                    if mask_u[y, x] == 1.:
                        u_right = data_u[y, x];

                    # v_down
                    v_down = 0
                    if mask_v[y - 1, x] == 1.:
                        v_down = data_v[y - 1, x];

                    # v_up
                    v_up = 0
                    if mask_v[y, x] == 1.:
                        v_up = data_v[y, x];

                    # 1.3 On calcule la demi-somme
                    u[y, x] = 0.5 * (u_left + u_right)
                    v[y, x] = 0.5 * (v_down + v_up)

                    # 1.4 On applique la rotation
                    u_rot[y, x] = u[y, x] * self.gridrotcos_t[y, x] + v[y, x] * self.gridrotsin_t[y, x]
                    v_rot[y, x] = -u[y, x] * self.gridrotsin_t[y, x] + v[y, x] * self.gridrotcos_t[y, x]

        # 2. On duplique les point sur les bords.
        # bottom
        u_rot[0, 0:self.xmax] = u_rot[1, 0:self.xmax]
        v_rot[0, 0:self.xmax] = v_rot[1, 0:self.xmax]
        # up
        u_rot[self.ymax - 1, 0:self.xmax] = u_rot[self.ymax - 2, 0:self.xmax]
        v_rot[self.ymax - 1, 0:self.xmax] = v_rot[self.ymax - 2, 0:self.xmax]

        # left
        u_rot[0:self.ymax, 0] = u_rot[0:self.ymax, 1]
        v_rot[0:self.ymax, 0] = v_rot[0:self.ymax, 1]
        # right
        u_rot[0:self.ymax, self.xmax - 1] = u_rot[0:self.ymax, self.xmax - 2]
        v_rot[0:self.ymax, self.xmax - 1] = v_rot[0:self.ymax, self.xmax - 2]

        return [u_rot, v_rot]

    def read_variable_waves_momentum_flux_to_ocean_at_time(self, index_t):
        index_z = self.zmax - 1
        mask_t = self.read_variable_2D_sea_binary_mask();
        mask_u = self.grid.variables["mask_u"][index_z][:];
        mask_v = self.grid.variables["mask_v"][index_z][:];
        data_u = self.ncfile.variables["twox"][index_t][::]
        data_v = self.ncfile.variables["twoy"][index_t][::]

        u = np.zeros([self.ymax, self.xmax])
        u[:] = np.NAN
        v = np.zeros([self.ymax, self.xmax])
        v[:] = np.NAN
        u_rot = np.zeros([self.ymax, self.xmax])
        u_rot[:] = np.NAN
        v_rot = np.zeros([self.ymax, self.xmax])
        v_rot[:] = np.NAN

        # 1. On calcule les point à l'intérieur du domaine en excluant les bords
        for y in range(1, self.ymax - 1):
            for x in range(1, self.xmax - 1):

                if mask_t[y, x] == 1.:

                    # 1.2 On récupère les valeurs aux point encadrant X pour faire la demi-somme
                    ##############################
                    #           v_up
                    #
                    #   u_left   X     u_right
                    #
                    #         v_bottom
                    #############################

                    # u_left
                    u_left = 0
                    if mask_u[y, x - 1] == 1.:
                        u_left = data_u[y, x - 1];

                    # u_right
                    u_right = 0
                    if mask_u[y, x] == 1.:
                        u_right = data_u[y, x];

                    # v_down
                    v_down = 0
                    if mask_v[y - 1, x] == 1.:
                        v_down = data_v[y - 1, x];

                    # v_up
                    v_up = 0
                    if mask_v[y, x] == 1.:
                        v_up = data_v[y, x];

                    # 1.3 On calcule la demi-somme
                    u[y, x] = 0.5 * (u_left + u_right)
                    v[y, x] = 0.5 * (v_down + v_up)

                    # 1.4 On applique la rotation
                    u_rot[y, x] = u[y, x] * self.gridrotcos_t[y, x] + v[y, x] * self.gridrotsin_t[y, x]
                    v_rot[y, x] = -u[y, x] * self.gridrotsin_t[y, x] + v[y, x] * self.gridrotcos_t[y, x]

        # 2. On duplique les point sur les bords.
        # bottom
        u_rot[0, 0:self.xmax] = u_rot[1, 0:self.xmax]
        v_rot[0, 0:self.xmax] = v_rot[1, 0:self.xmax]
        # up
        u_rot[self.ymax - 1, 0:self.xmax] = u_rot[self.ymax - 2, 0:self.xmax]
        v_rot[self.ymax - 1, 0:self.xmax] = v_rot[self.ymax - 2, 0:self.xmax]

        # left
        u_rot[0:self.ymax, 0] = u_rot[0:self.ymax, 1]
        v_rot[0:self.ymax, 0] = v_rot[0:self.ymax, 1]
        # right
        u_rot[0:self.ymax, self.xmax - 1] = u_rot[0:self.ymax, self.xmax - 2]
        v_rot[0:self.ymax, self.xmax - 1] = v_rot[0:self.ymax, self.xmax - 2]
        return [u_rot, v_rot]

    #################
    # METEO
    # 2D
    #################


    #################
    # METEO
    # Sea surface
    #################

    def read_variable_wind_stress_at_time(self, index_t):
        index_z = self.zmax - 1
        mask_t = self.read_variable_2D_sea_binary_mask();
        mask_u = self.grid.variables["mask_u"][index_z][:];
        mask_v = self.grid.variables["mask_v"][index_z][:];
        data_u = self.ncfile.variables["wstress_u"][index_t][::]
        data_v = self.ncfile.variables["wstress_v"][index_t][::]

        u = np.zeros([self.ymax, self.xmax])
        u[:] = np.NAN
        v = np.zeros([self.ymax, self.xmax])
        v[:] = np.NAN
        u_rot = np.zeros([self.ymax, self.xmax])
        u_rot[:] = np.NAN
        v_rot = np.zeros([self.ymax, self.xmax])
        v_rot[:] = np.NAN

        # 1. On calcule les point à l'intérieur du domaine en excluant les bords
        for y in range(1, self.ymax - 1):
            for x in range(1, self.xmax - 1):

                if mask_t[y, x] == 1.:

                    # 1.2 On récupère les valeurs aux point encadrant X pour faire la demi-somme
                    ##############################
                    #           v_up
                    #
                    #   u_left   X     u_right
                    #
                    #         v_bottom
                    #############################

                    # u_left
                    u_left = 0
                    if mask_u[y, x - 1] == 1.:
                        u_left = data_u[y, x - 1];

                    # u_right
                    u_right = 0
                    if mask_u[y, x] == 1.:
                        u_right = data_u[y, x];

                    # v_down
                    v_down = 0
                    if mask_v[y - 1, x] == 1.:
                        v_down = data_v[y - 1, x];

                    # v_up
                    v_up = 0
                    if mask_v[y, x] == 1.:
                        v_up = data_v[y, x];

                    # 1.3 On calcule la demi-somme
                    u[y, x] = 0.5 * (u_left + u_right)
                    v[y, x] = 0.5 * (v_down + v_up)

                    # 1.4 On applique la rotation
                    u_rot[y, x] = u[y, x] * self.gridrotcos_t[y, x] + v[y, x] * self.gridrotsin_t[y, x]
                    v_rot[y, x] = -u[y, x] * self.gridrotsin_t[y, x] + v[y, x] * self.gridrotcos_t[y, x]

        # 2. On duplique les point sur les bords.
        # bottom
        u_rot[0, 0:self.xmax] = u_rot[1, 0:self.xmax]
        v_rot[0, 0:self.xmax] = v_rot[1, 0:self.xmax]
        # up
        u_rot[self.ymax - 1, 0:self.xmax] = u_rot[self.ymax - 2, 0:self.xmax]
        v_rot[self.ymax - 1, 0:self.xmax] = v_rot[self.ymax - 2, 0:self.xmax]

        # left
        u_rot[0:self.ymax, 0] = u_rot[0:self.ymax, 1]
        v_rot[0:self.ymax, 0] = v_rot[0:self.ymax, 1]
        # right
        u_rot[0:self.ymax, self.xmax - 1] = u_rot[0:self.ymax, self.xmax - 2]
        v_rot[0:self.ymax, self.xmax - 1] = v_rot[0:self.ymax, self.xmax - 2]

        return [u_rot, v_rot]

    #################
    # METEO
    # At 10 m
    #################

    def read_variable_wind_10m_at_time(self, index_t):
        mask_t = self.read_variable_2D_sea_binary_mask();
        lon_t = self.read_axis_x();
        lat_t = self.read_axis_y();
        data_u = self.ncfile.variables["uwind_t"][index_t][::]
        data_v = self.ncfile.variables["vwind_t"][index_t][::]

        u_rot = np.zeros([self.ymax, self.xmax])
        u_rot[:] = np.NAN
        v_rot = np.zeros([self.ymax, self.xmax])
        v_rot[:] = np.NAN

        # 1. On calcule les point à l'intérieur du domaine en excluant les bords
        for y in range(1, self.ymax - 1):
            for x in range(1, self.xmax - 1):

                if mask_t[y, x] == 1.:
                    # 1.4 On applique la rotation
                    u_rot[y, x] = data_u[y, x] * self.gridrotcos_t[y, x] + data_v[y, x] * self.gridrotsin_t[y, x]
                    v_rot[y, x] = -data_u[y, x] * self.gridrotsin_t[y, x] + data_v[y, x] * self.gridrotcos_t[y, x]

        # 2. On duplique les point sur les bords.
        # bottom
        u_rot[0, 0:self.xmax] = u_rot[1, 0:self.xmax]
        v_rot[0, 0:self.xmax] = v_rot[1, 0:self.xmax]
        # up
        u_rot[self.ymax - 1, 0:self.xmax] = u_rot[self.ymax - 2, 0:self.xmax]
        v_rot[self.ymax - 1, 0:self.xmax] = v_rot[self.ymax - 2, 0:self.xmax]

        # left
        u_rot[0:self.ymax, 0] = u_rot[0:self.ymax, 1]
        v_rot[0:self.ymax, 0] = v_rot[0:self.ymax, 1]
        # right
        u_rot[0:self.ymax, self.xmax - 1] = u_rot[0:self.ymax, self.xmax - 2]
        v_rot[0:self.ymax, self.xmax - 1] = v_rot[0:self.ymax, self.xmax - 2]

        return [u_rot, v_rot]
        
    # Others
    def read_variable_Ha(self):
        return self.ncfile.variables["Ha"][:]

    def read_variable_bathy_ssh_at_time(self,index_t):
        return self.ncfile.variables["hssh"][index_t][:]


