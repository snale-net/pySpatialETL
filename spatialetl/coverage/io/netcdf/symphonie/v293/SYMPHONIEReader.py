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

import numpy as np

from spatialetl.coverage.io.netcdf.symphonie.SYMPHONIEReader import SYMPHONIEReader as AbstractSYMPHONIEReader
from spatialetl.exception.VariableNameError import VariableNameError
from spatialetl.utils.VariableDefinition import VariableDefinition
from spatialetl.utils.logger import logging


class SYMPHONIEReader(AbstractSYMPHONIEReader):
    """
La classe SymphonieReader permet de lire les données du format Symphonie

@param  myGrid : lien vers le fichier de grille (que l'on trouve dans le RDIR/tmp/grid.nc)
@param myFile : lien vers le fichier de données (que l'on trouve dans GRAPHIQUES)
"""

    def __init__(self,myGrid, myFile):
        AbstractSYMPHONIEReader.__init__(self, myGrid, myFile);

    def read_variable_2D_sea_binary_mask(self,xmin,xmax,ymin,ymax):
        index_z = self.get_z_size() - 1 # At surface level
        try:
            if "mask_t" in self.grid.variables:
                return np.ma.filled(self.grid.variables["mask_t"][ymin:ymax, xmin:xmax], fill_value=np.nan)
        except Exception as ex:
            logging.debug("Error '" + str(ex) + "'")
            raise (VariableNameError("SymphonieReader", "An error occured : '" + str(ex) + "'", 1000))

        logging.debug("No variables found for '" + str(VariableDefinition.LONG_NAME['2d_sea_binary_mask']) + "'")
        raise (VariableNameError("SymphonieReader",
                                 "No variables found for '" + str(VariableDefinition.LONG_NAME['2d_sea_binary_mask']) + "'",
                                 1000))

    def read_variable_3D_sea_binary_mask(self,xmin,xmax,ymin,ymax):
        try:
            if "mask_t" in self.grid.variables:
                return np.ma.filled(self.grid.variables["mask_t"][:,ymin:ymax, xmin:xmax], fill_value=np.nan)
        except Exception as ex:
            logging.debug("Error '" + str(ex) + "'")
            raise (VariableNameError("SymphonieReader", "An error occured : '" + str(ex) + "'", 1000))

        logging.debug("No variables found for '" + str(VariableDefinition.LONG_NAME['2d_sea_binary_mask']) + "'")
        raise (VariableNameError("SymphonieReader",
                                 "No variables found for '" + str(
                                     VariableDefinition.LONG_NAME['2d_sea_binary_mask']) + "'",
                                 1000))

    #################
    # HYDRO
    # Sea Surface
    #################

    def read_variable_sea_water_velocity_at_sea_water_surface_at_time(self, index_t, xmin, xmax, ymin, ymax):
        self.open_file(index_t)
        index_z = self.get_z_size() - 1
        xmin_overlap, xmax_overlap, ymin_overlap, ymax_overlap, new_xmin, new_xmax, new_ymin, new_ymax = self.compute_overlap_indexes(
            xmin, xmax, ymin, ymax)

        mask_t = self.grid.variables["mask_t"][ymin_overlap:ymax_overlap, xmin_overlap:xmax_overlap];
        mask_u = self.grid.variables["mask_u"][ymin_overlap:ymax_overlap, xmin_overlap:xmax_overlap];
        mask_v = self.grid.variables["mask_v"][ymin_overlap:ymax_overlap, xmin_overlap:xmax_overlap];

        if self.gridrotcos_t is None and self.gridrotsin_t is None:
            self.compute_rot()

        rotcos = self.gridrotcos_t[ymin_overlap:ymax_overlap, xmin_overlap:xmax_overlap];
        rotsin = self.gridrotsin_t[ymin_overlap:ymax_overlap, xmin_overlap:xmax_overlap];

        if "vel_u" in self.ncfile.variables:
            data_u = np.ma.filled(
                self.ncfile.variables["vel_u"][0, index_z, ymin_overlap:ymax_overlap, xmin_overlap:xmax_overlap],
                fill_value=np.nan)
        if "u" in self.ncfile.variables:
            data_u = np.ma.filled(
                self.ncfile.variables["u"][0, index_z, ymin_overlap:ymax_overlap, xmin_overlap:xmax_overlap],
                fill_value=np.nan)
        if "vel_v" in self.ncfile.variables:
            data_v = np.ma.filled(
                self.ncfile.variables["vel_v"][0, index_z, ymin_overlap:ymax_overlap, xmin_overlap:xmax_overlap],
                fill_value=np.nan)
        if "v" in self.ncfile.variables:
            data_v = np.ma.filled(
                self.ncfile.variables["v"][0, index_z, ymin_overlap:ymax_overlap, xmin_overlap:xmax_overlap],
                fill_value=np.nan)

        u_rot, v_rot = self.compute_vector_rotation(data_u, data_v, rotcos, rotsin, mask_t, mask_u, mask_v)

        if "wetmask_t" in self.ncfile.variables:
            u_rot[self.ncfile.variables["wetmask_t"][0, ymin_overlap:ymax_overlap,
                  xmin_overlap:xmax_overlap] == 0] = np.nan
            v_rot[self.ncfile.variables["wetmask_t"][0, ymin_overlap:ymax_overlap,
                  xmin_overlap:xmax_overlap] == 0] = np.nan

        return [u_rot[new_ymin:new_ymax, new_xmin:new_xmax], v_rot[new_ymin:new_ymax, new_xmin:new_xmax]]

    #################
    # HYDRO
    # Ground level
    #################

    #################
    # HYDRO
    # 2D
    #################

    def read_variable_barotropic_sea_water_velocity_at_time(self, index_t, xmin, xmax, ymin, ymax):
        try:
            self.open_file(index_t)
            index_z = self.get_z_size() - 1
            xmin_overlap, xmax_overlap, ymin_overlap, ymax_overlap, new_xmin, new_xmax, new_ymin, new_ymax = self.compute_overlap_indexes(
                xmin, xmax, ymin, ymax)

            mask_t = self.grid.variables["mask_t"][ymin_overlap:ymax_overlap, xmin_overlap:xmax_overlap];
            mask_u = self.grid.variables["mask_u"][ymin_overlap:ymax_overlap, xmin_overlap:xmax_overlap];
            mask_v = self.grid.variables["mask_v"][ymin_overlap:ymax_overlap, xmin_overlap:xmax_overlap];

            if self.gridrotcos_t is None and self.gridrotsin_t is None:
                self.compute_rot()

            rotcos = self.gridrotcos_t[ymin_overlap:ymax_overlap, xmin_overlap:xmax_overlap];
            rotsin = self.gridrotsin_t[ymin_overlap:ymax_overlap, xmin_overlap:xmax_overlap];

            if "velbar_u" in self.ncfile.variables:
                data_u = np.ma.filled(
                    self.ncfile.variables["velbar_u"][0, ymin_overlap:ymax_overlap, xmin_overlap:xmax_overlap],
                    fill_value=np.nan)
            if "velbar_v" in self.ncfile.variables:
                data_v = np.ma.filled(
                    self.ncfile.variables["velbar_v"][0, ymin_overlap:ymax_overlap, xmin_overlap:xmax_overlap],
                    fill_value=np.nan)
            u_rot, v_rot = self.compute_vector_rotation(data_u, data_v, rotcos, rotsin, mask_t, mask_u, mask_v)

            if "wetmask_t" in self.ncfile.variables:
                u_rot[self.ncfile.variables["wetmask_t"][0, ymin_overlap:ymax_overlap,
                      xmin_overlap:xmax_overlap] == 0] = np.nan
                v_rot[self.ncfile.variables["wetmask_t"][0, ymin_overlap:ymax_overlap,
                      xmin_overlap:xmax_overlap] == 0] = np.nan

            return [u_rot[new_ymin:new_ymax, new_xmin:new_xmax], v_rot[new_ymin:new_ymax, new_xmin:new_xmax]]

        except Exception as ex:
            logging.debug("Error '" + str(ex) + "'")
            raise (VariableNameError("SymphonieReader", "An error occured : '" + str(ex) + "'", 1000))

        logging.debug("No variables found for \'Barotropic Sea Water Velocity\'")
        raise (VariableNameError("SymphonieReader",
                                 "No variables found for \'Barotropic Sea Water Velocity\'",
                                 1000))

    #################
    # HYDRO
    # 3D
    #################

    def read_variable_baroclinic_sea_water_velocity_at_time_and_depth(self, index_t, index_z, xmin, xmax, ymin, ymax):
        try:
            self.open_file(index_t)
            xmin_overlap, xmax_overlap, ymin_overlap, ymax_overlap, new_xmin, new_xmax, new_ymin, new_ymax = self.compute_overlap_indexes(
                xmin, xmax, ymin, ymax)

            mask_t = self.grid.variables["mask_t"][ymin_overlap:ymax_overlap, xmin_overlap:xmax_overlap];
            mask_u = self.grid.variables["mask_u"][ymin_overlap:ymax_overlap, xmin_overlap:xmax_overlap];
            mask_v = self.grid.variables["mask_v"][ymin_overlap:ymax_overlap, xmin_overlap:xmax_overlap];

            if self.gridrotcos_t is None and self.gridrotsin_t is None:
                self.compute_rot()

            rotcos = self.gridrotcos_t[ymin_overlap:ymax_overlap, xmin_overlap:xmax_overlap];
            rotsin = self.gridrotsin_t[ymin_overlap:ymax_overlap, xmin_overlap:xmax_overlap];

            if "vel_u" in self.ncfile.variables:
                data_u = np.ma.filled(
                    self.ncfile.variables["vel_u"][0, index_z, ymin_overlap:ymax_overlap,
                    xmin_overlap:xmax_overlap],
                    fill_value=np.nan)
            if "vel_v" in self.ncfile.variables:
                data_v = np.ma.filled(
                    self.ncfile.variables["vel_v"][0, index_z, ymin_overlap:ymax_overlap,
                    xmin_overlap:xmax_overlap],
                    fill_value=np.nan)

            u_rot, v_rot = self.compute_vector_rotation(data_u, data_v, rotcos, rotsin, mask_t, mask_u, mask_v)

            if "wetmask_t" in self.ncfile.variables:
                u_rot[self.ncfile.variables["wetmask_t"][0, ymin_overlap:ymax_overlap,
                      xmin_overlap:xmax_overlap] == 0] = np.nan
                v_rot[self.ncfile.variables["wetmask_t"][0, ymin_overlap:ymax_overlap,
                      xmin_overlap:xmax_overlap] == 0] = np.nan

            return [u_rot[new_ymin:new_ymax, new_xmin:new_xmax], v_rot[new_ymin:new_ymax, new_xmin:new_xmax]]

        except Exception as ex:
            logging.debug("Error '" + str(ex) + "'")
            raise (VariableNameError("SymphonieReader", "An error occured : '" + str(ex) + "'", 1000))

        logging.debug("No variables found for \'Baroclinic Sea Water Velocity\'")
        raise (VariableNameError("SymphonieReader",
                                 "No variables found for \'Baroclinic Sea Water Velocity\'",
                                 1000))

    #################
    # WAVES
    # Sea Surface
    #################

    #################
    # WAVES
    # Momentum flux
    #################

    #################
    # METEO
    # 2D
    #################

    #################
    # METEO
    # Sea surface
    #################
    def read_variable_wind_stress_at_time(self, index_t, xmin, xmax, ymin, ymax):
        try:
            self.open_file(index_t)
            index_z = self.get_z_size() - 1
            xmin_overlap, xmax_overlap, ymin_overlap, ymax_overlap, new_xmin, new_xmax, new_ymin, new_ymax = self.compute_overlap_indexes(
                xmin, xmax, ymin, ymax)

            mask_t = self.grid.variables["mask_t"][ymin_overlap:ymax_overlap, xmin_overlap:xmax_overlap];
            mask_u = self.grid.variables["mask_u"][ymin_overlap:ymax_overlap, xmin_overlap:xmax_overlap];
            mask_v = self.grid.variables["mask_v"][ymin_overlap:ymax_overlap, xmin_overlap:xmax_overlap];

            if self.gridrotcos_t is None and self.gridrotsin_t is None:
                self.compute_rot()

            rotcos = self.gridrotcos_t[ymin_overlap:ymax_overlap, xmin_overlap:xmax_overlap];
            rotsin = self.gridrotsin_t[ymin_overlap:ymax_overlap, xmin_overlap:xmax_overlap];

            if "wstress_u" in self.ncfile.variables:
                data_u = np.ma.filled(
                    self.ncfile.variables["wstress_u"][0, ymin_overlap:ymax_overlap, xmin_overlap:xmax_overlap],
                    fill_value=np.nan)
            if "wstress_v" in self.ncfile.variables:
                data_v = np.ma.filled(
                    self.ncfile.variables["wstress_v"][0, ymin_overlap:ymax_overlap, xmin_overlap:xmax_overlap],
                    fill_value=np.nan)

            u_rot, v_rot = self.compute_vector_rotation(data_u, data_v, rotcos, rotsin, mask_t, mask_u, mask_v)

            if "wetmask_t" in self.ncfile.variables:
                u_rot[self.ncfile.variables["wetmask_t"][0, ymin_overlap:ymax_overlap,
                      xmin_overlap:xmax_overlap] == 0] = np.nan
                v_rot[self.ncfile.variables["wetmask_t"][0, ymin_overlap:ymax_overlap,
                      xmin_overlap:xmax_overlap] == 0] = np.nan

            return [u_rot[new_ymin:new_ymax, new_xmin:new_xmax], v_rot[new_ymin:new_ymax, new_xmin:new_xmax]]

        except Exception as ex:
            logging.debug("Error '" + str(ex) + "'")
            raise (VariableNameError("SymphonieReader", "An error occured : '" + str(ex) + "'", 1000))

        logging.debug("No variables found for \'Wind Stress\'")
        raise (VariableNameError("SymphonieReader",
                                 "No variables found for \'Wind Stress\'",
                                 1000))

        #################
        # METEO
        # At 10 m
        #################

    #################
    # METEO
    # At 10 m
    #################
    def read_variable_wind_10m_at_time(self, index_t, xmin, xmax, ymin, ymax):
        try:
            self.open_file(index_t)
            index_z = self.get_z_size() - 1
            xmin_overlap, xmax_overlap, ymin_overlap, ymax_overlap, new_xmin, new_xmax, new_ymin, new_ymax = self.compute_overlap_indexes(
                xmin, xmax, ymin, ymax)

            mask_t = self.grid.variables["mask_t"][ymin_overlap:ymax_overlap, xmin_overlap:xmax_overlap];
            mask_u = self.grid.variables["mask_u"][ymin_overlap:ymax_overlap, xmin_overlap:xmax_overlap];
            mask_v = self.grid.variables["mask_v"][ymin_overlap:ymax_overlap, xmin_overlap:xmax_overlap];

            if self.gridrotcos_t is None and self.gridrotsin_t is None:
                self.compute_rot()

            rotcos = self.gridrotcos_t[ymin_overlap:ymax_overlap, xmin_overlap:xmax_overlap];
            rotsin = self.gridrotsin_t[ymin_overlap:ymax_overlap, xmin_overlap:xmax_overlap];

            if "uwind_t" in self.ncfile.variables:
                data_u = np.ma.filled(
                    self.ncfile.variables["uwind_t"][0, ymin_overlap:ymax_overlap, xmin_overlap:xmax_overlap],
                    fill_value=np.nan)
            if "vwind_t" in self.ncfile.variables:
                data_v = np.ma.filled(
                    self.ncfile.variables["vwind_t"][0, ymin_overlap:ymax_overlap, xmin_overlap:xmax_overlap],
                    fill_value=np.nan)

            u_rot, v_rot = self.compute_vector_rotation(data_u, data_v, rotcos, rotsin, mask_t, mask_u, mask_v)

            if "wetmask_t" in self.ncfile.variables:
                u_rot[self.ncfile.variables["wetmask_t"][0, ymin_overlap:ymax_overlap,
                      xmin_overlap:xmax_overlap] == 0] = np.nan
                v_rot[self.ncfile.variables["wetmask_t"][0, ymin_overlap:ymax_overlap,
                      xmin_overlap:xmax_overlap] == 0] = np.nan

            return [u_rot[new_ymin:new_ymax, new_xmin:new_xmax], v_rot[new_ymin:new_ymax, new_xmin:new_xmax]]

        except Exception as ex:
            logging.debug("Error '" + str(ex) + "'")
            raise (VariableNameError("SymphonieReader", "An error occured : '" + str(ex) + "'", 1000))

        logging.debug("No variables found for \'Wind 10m\'")
        raise (VariableNameError("SymphonieReader",
                                 "No variables found for \'Wind 10m\'",
                                 1000))