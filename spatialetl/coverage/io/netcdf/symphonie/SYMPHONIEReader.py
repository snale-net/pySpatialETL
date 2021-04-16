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

import glob
import os
import re

import cftime
import numpy as np
from netCDF4 import Dataset, num2date

from spatialetl.coverage.TimeCoverage import TimeCoverage
from spatialetl.coverage.io.CoverageReader import CoverageReader
from spatialetl.exception.VariableNameError import VariableNameError
from spatialetl.utils.VariableDefinition import VariableDefinition
from spatialetl.utils.logger import logging
from spatialetl.utils.path import path_leaf


class SYMPHONIEReader(CoverageReader):
    """
La classe SymphonieReader permet de lire les données du format Symphonie

@param  myGrid : lien vers le fichier de grille (que l'on trouve dans le RDIR/tmp/grid.nc)
@param myFile : lien vers le fichier de données (que l'on trouve dans GRAPHIQUES)
"""
    HORIZONTAL_OVERLAPING_SIZE = 2

    def __init__(self,myGrid, myFile=None):
        CoverageReader.__init__(self,myGrid);

        self.grid = Dataset(self.filename, 'r')
        self.gridrotcos_t = None
        self.gridrotsin_t = None

        if myFile is not None:
            if os.path.isfile(myFile):
                self.files= [myFile]
            elif os.path.isdir(myFile):
                self.files = sorted(glob.glob(os.path.join(myFile, "*.nc")))
            elif myFile.endswith("*"):
                self.files = sorted(glob.glob(myFile + ".nc"))
            else:
                raise ValueError("Unable to decode file " + str(myFile))

            self.ncfile = Dataset(self.files[0], 'r')
            self.last_opened_t_index = 0
        else:
            self.files = []
            self.ncfile = Dataset(self.filename, 'r')
            self.last_opened_t_index = 0

        self.t_size = len(self.files)
        self.times = []

        for file in self.files:
            groups = re.search("^([0-9]{4})([0-9]{2})([0-9]{2})\_([0-9]{2})([0-9]{2})([0-9]{2}).*.nc$",
                               path_leaf(file))
            if groups:
                current_time = cftime.datetime(int(groups.group(1)), int(groups.group(2)), int(groups.group(3)),
                                               int(groups.group(4)), int(groups.group(5)),
                                               int(groups.group(6)))
                self.times.append(current_time)
            else:
                if os.path.isfile(file):
                    try:
                        current_file = Dataset(file, 'r')

                        if "time" in current_file.variables and np.shape(current_file.variables["time"]) == (1,):
                            current_time = num2date(current_file.variables["time"][0],
                                     units=current_file.variables["time"].units.replace('from', 'since').replace('jan',
                                                                                                                    '01').replace(
                                         'feb', '02').replace('mar', '03').replace('apr', '04').replace('may',
                                                                                                        '05').replace('jun',
                                                                                                                      '06').replace(
                                         'jul', '07').replace('aug', '08').replace('sep', '09').replace('oct',
                                                                                                        '10').replace('nov',
                                                                                                                      '11').replace(
                                         'dec', '12'), calendar=current_file.variables['time'].calendar)
                            self.times.append(current_time)
                        else:
                            raise ValueError("No variable time found or multiple time records found in the same file")
                    except Exception as ex:
                        raise ValueError("Unable to decode time records in file "+str(file)+ ":"+str(ex))

        if len(self.times) == 0:
            raise ValueError("No time records found")

    def open_file(self, index_t):
        if index_t != self.last_opened_t_index:
            self.close()
            self.ncfile = Dataset(self.files[index_t])

    def close(self):
        self.ncfile.close()

    def compute_rot(self):

        logging.debug("[SymphonieReader] Compute grid rotation matrix...")

        lon_t = np.ma.filled(self.grid.variables['longitude_t'][:], fill_value=np.nan)
        lat_t = np.ma.filled(self.grid.variables['latitude_t'][:], fill_value=np.nan)

        self.gridrotcos_t = np.zeros([self.get_y_size(), self.get_x_size()])
        self.gridrotsin_t = np.zeros([self.get_y_size(), self.get_x_size()])

        for y in range(1, self.get_y_size() - 1):
            for x in range(1, self.get_x_size() - 1):

                x1 = (lon_t[y, x + 1] - lon_t[y, x - 1]) * np.pi / 180.
                if (x1 < -np.pi): x1 = x1 + 2. * np.pi
                if (x1 > np.pi): x1 = x1 - 2. * np.pi
                x0 = -np.arctan2((lat_t[y, x + 1] - lat_t[y, x - 1]) * np.pi / 180.,
                                 x1 * np.cos(lat_t[y, x] * np.pi / 180.))
                self.gridrotcos_t[y, x] = np.cos(x0)
                self.gridrotsin_t[y, x] = np.sin(x0)

    def compute_to_tracer(self, data_u, data_v, mask_t, mask_u, mask_v):

        x_size = np.shape(mask_t)[1]
        y_size = np.shape(mask_t)[0]

        u_t = np.zeros([y_size, x_size])
        u_t[:] = np.nan
        v_t = np.zeros([y_size, x_size])
        v_t[:] = np.nan

        # 1. On calcule les point à l'intérieur du domaine en excluant les bords
        for y in range(1, y_size - 1):
            for x in range(1, x_size - 1):

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
                    u_t[y, x] = 0.5 * (u_left + u_right)
                    v_t[y, x] = 0.5 * (v_down + v_up)

        # 2. On duplique les point sur les bords.
        # bottom
        u_t[0, 0:x_size] = u_t[1, 0:x_size]
        v_t[0, 0:x_size] = v_t[1, 0:x_size]
        # up
        u_t[y_size - 1, 0:x_size] = u_t[y_size - 2, 0:x_size]
        v_t[y_size - 1, 0:x_size] = v_t[y_size - 2, 0:x_size]

        # left
        u_t[0:y_size, 0] = u_t[0:y_size, 1]
        v_t[0:y_size, 0] = v_t[0:y_size, 1]
        # right
        u_t[0:y_size, x_size - 1] = u_t[0:y_size, x_size - 2]
        v_t[0:y_size, x_size - 1] = v_t[0:y_size, x_size - 2]

        return u_t, v_t

    def compute_vector_rotation(self, data_u, data_v, rotcos, rotsin, mask_t, mask_u, mask_v):

        x_size = np.shape(mask_t)[1]
        y_size = np.shape(mask_t)[0]

        u = np.zeros([y_size, x_size])
        u[:] = np.nan
        v = np.zeros([y_size, x_size])
        v[:] = np.nan
        u_rot = np.zeros([y_size, x_size])
        u_rot[:] = np.nan
        v_rot = np.zeros([y_size, x_size])
        v_rot[:] = np.nan

        # 1. On calcule les point à l'intérieur du domaine en excluant les bords
        for y in range(1, y_size - 1):
            for x in range(1, x_size - 1):

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
                    u_rot[y, x] = u[y, x] * rotcos[y, x] + v[y, x] * rotsin[y, x]
                    v_rot[y, x] = -u[y, x] * rotsin[y, x] + v[y, x] * rotcos[y, x]

        # 2. On duplique les point sur les bords.
        # bottom
        u_rot[0, 0:x_size] = u_rot[1, 0:x_size]
        v_rot[0, 0:x_size] = v_rot[1, 0:x_size]
        # up
        u_rot[y_size - 1, 0:x_size] = u_rot[y_size - 2, 0:x_size]
        v_rot[y_size - 1, 0:x_size] = v_rot[y_size - 2, 0:x_size]

        # left
        u_rot[0:y_size, 0] = u_rot[0:y_size, 1]
        v_rot[0:y_size, 0] = v_rot[0:y_size, 1]
        # right
        u_rot[0:y_size, x_size - 1] = u_rot[0:y_size, x_size - 2]
        v_rot[0:y_size, x_size - 1] = v_rot[0:y_size, x_size - 2]

        return u_rot, v_rot

    def compute_overlap_indexes(self, xmin, xmax, ymin, ymax):

        xmin_overlap = max(0, xmin - SYMPHONIEReader.HORIZONTAL_OVERLAPING_SIZE)
        new_xmin = SYMPHONIEReader.HORIZONTAL_OVERLAPING_SIZE
        if xmin_overlap == 0:
            new_xmin = xmin

        xmax_overlap = min(self.get_x_size(), xmax + SYMPHONIEReader.HORIZONTAL_OVERLAPING_SIZE)
        new_xmax = xmax_overlap - xmin_overlap - SYMPHONIEReader.HORIZONTAL_OVERLAPING_SIZE
        if xmax_overlap == self.get_x_size():
            new_xmax = self.get_x_size()

        ymin_overlap = max(0, ymin - SYMPHONIEReader.HORIZONTAL_OVERLAPING_SIZE)
        new_ymin = SYMPHONIEReader.HORIZONTAL_OVERLAPING_SIZE
        if ymin_overlap == 0:
            new_ymin = ymin

        ymax_overlap = min(self.get_y_size(), ymax + SYMPHONIEReader.HORIZONTAL_OVERLAPING_SIZE)
        new_ymax = ymax_overlap - ymin_overlap - SYMPHONIEReader.HORIZONTAL_OVERLAPING_SIZE
        if ymax_overlap == self.get_y_size():
            new_ymax = self.get_y_size()

        return xmin_overlap, xmax_overlap, ymin_overlap, ymax_overlap, new_xmin, new_xmax, new_ymin, new_ymax

    def is_regular_grid(self):
        return False

    def get_x_size(self):
        try:
            if "longitude_t" in self.grid.variables:
                return np.shape(self.grid.variables['longitude_t'][:])[1];
            else:
                logging.debug(
                    "No variables found for '" + str(VariableDefinition.LONG_NAME['longitude']) + "'")
                raise (VariableNameError("SymphonieReader",
                                         "No variables found for '" + str(
                                             VariableDefinition.LONG_NAME['longitude']) + "'",
                                         1000))
        except Exception as ex:
            logging.debug("Error '" + str(ex) + "'")
            raise (VariableNameError("SymphonieReader", "An error occured : '" + str(ex) + "'", 1000))

    def get_y_size(self):
        try:
            if "latitude_t" in self.grid.variables:
                return np.shape(self.grid.variables['latitude_t'][:])[0];
            else:
                logging.debug(
                    "No variables found for '" + str(VariableDefinition.LONG_NAME['latitude']) + "'")
                raise (VariableNameError("SymphonieReader",
                                         "No variables found for '" + str(
                                             VariableDefinition.LONG_NAME['latitude']) + "'",
                                         1000))
        except Exception as ex:
            logging.debug("Error '" + str(ex) + "'")
            raise (VariableNameError("SymphonieReader", "An error occured : '" + str(ex) + "'", 1000))

    def get_z_size(self):
        try:
            if "depth_t" in self.grid.variables:
                return np.shape(self.grid.variables['depth_t'][:])[0];
            else:
                logging.debug(
                    "No variables found for '" + str(VariableDefinition.LONG_NAME['depth']) + "'")
                raise (VariableNameError("SymphonieReader",
                                         "No variables found for '" + str(
                                             VariableDefinition.LONG_NAME['depth']) + "'",
                                         1000))
        except Exception as ex:
            logging.debug("Error '" + str(ex) + "'")
            raise (VariableNameError("SymphonieReader", "An error occured : '" + str(ex) + "'", 1000))

    def get_t_size(self):
        return self.t_size

    def read_axis_x(self, xmin, xmax, ymin, ymax):
        try:
            if "longitude_t" in self.grid.variables:
                return self.grid.variables['longitude_t'][ymin:ymax, xmin:xmax]
            else:
                logging.debug(
                    "No variables found for '" + str(VariableDefinition.LONG_NAME['longitude']) + "'")
                raise (VariableNameError("SymphonieReader",
                                         "No variables found for '" + str(
                                             VariableDefinition.LONG_NAME['longitude']) + "'",
                                         1000))
        except Exception as ex:
            logging.debug("Error '" + str(ex) + "'")
            raise (VariableNameError("SymphonieReader", "An error occured : '" + str(ex) + "'", 1000))

    def read_axis_y(self, xmin, xmax, ymin, ymax):
        try:
            if "latitude_t" in self.grid.variables:
                return self.grid.variables['latitude_t'][ymin:ymax, xmin:xmax]
            else:
                logging.debug(
                    "No variables found for '" + str(VariableDefinition.LONG_NAME['latitude']) + "'")
                raise (VariableNameError("SymphonieReader",
                                         "No variables found for '" + str(
                                             VariableDefinition.LONG_NAME['latitude']) + "'",
                                         1000))
        except Exception as ex:
            logging.debug("Error '" + str(ex) + "'")
            raise (VariableNameError("SymphonieReader", "An error occured : '" + str(ex) + "'", 1000))

    def read_axis_z(self, ):
        try:
            if "depth_t" in self.grid.variables:
                lev = self.grid.variables["depth_t"][::]
                lev[::] *= -1.0  # inverse la profondeur
                return lev
            else:
                logging.debug(
                    "No variables found for '" + str(VariableDefinition.LONG_NAME['depth']) + "'")
                raise (VariableNameError("SymphonieReader",
                                         "No variables found for '" + str(
                                             VariableDefinition.LONG_NAME['depth']) + "'",
                                         1000))
        except Exception as ex:
            logging.debug("Error '" + str(ex) + "'")
            raise (VariableNameError("SymphonieReader", "An error occured : '" + str(ex) + "'", 1000))

    def read_axis_t(self, tmin, tmax, timestamp):

        if timestamp == 1:
            return [(t - TimeCoverage.TIME_DATUM).total_seconds() \
                    for t in self.times[tmin:tmax]];
        else:
            return self.times[tmin:tmax]

    # Variables
    def read_variable_time(self,tmin, tmax, timestamp):
        return self.read_axis_t(tmin,tmax, timestamp=timestamp)

    def read_variable_longitude(self, xmin, xmax, ymin, ymax):
        return self.read_axis_x(xmin,xmax,ymin,ymax)

    def read_variable_latitude(self, xmin, xmax, ymin, ymax):
        return self.read_axis_y(xmin,xmax,ymin,ymax)

    def read_variable_2D_sea_binary_mask(self, xmin, xmax, ymin, ymax):
        index_z = self.get_z_size() - 1  # At surface level
        try:
            if "mask_t" in self.grid.variables:
                return np.ma.filled(self.grid.variables["mask_t"][index_z, ymin:ymax, xmin:xmax], fill_value=np.nan)
            else:
                logging.debug(
                    "No variables found for '" + str(VariableDefinition.LONG_NAME['2d_sea_binary_mask']) + "'")
                raise (VariableNameError("SymphonieReader",
                                         "No variables found for '" + str(
                                             VariableDefinition.LONG_NAME['2d_sea_binary_mask']) + "'",
                                         1000))
        except Exception as ex:
            logging.debug("Error '" + str(ex) + "'")
            raise (VariableNameError("SymphonieReader", "An error occured : '" + str(ex) + "'", 1000))

    def read_variable_3D_sea_binary_mask(self, xmin, xmax, ymin, ymax):
        try:
            if "mask_t" in self.grid.variables:
                return np.ma.filled(self.grid.variables["mask_t"][:, ymin:ymax, xmin:xmax], fill_value=np.nan)
            else:
                logging.debug(
                    "No variables found for '" + str(VariableDefinition.LONG_NAME['2d_sea_binary_mask']) + "'")
                raise (VariableNameError("SymphonieReader",
                                         "No variables found for '" + str(
                                             VariableDefinition.LONG_NAME['2d_sea_binary_mask']) + "'",
                                         1000))
        except Exception as ex:
            logging.debug("Error '" + str(ex) + "'")
            raise (VariableNameError("SymphonieReader", "An error occured : '" + str(ex) + "'", 1000))

    def read_variable_2D_wet_binary_mask_at_time(self, index_t, xmin, xmax, ymin, ymax):
        try:
            self.open_file(index_t)
            if "wetmask_t" in self.ncfile.variables:
                return np.ma.filled(self.ncfile.variables["wetmask_t"][0, ymin:ymax, xmin:xmax],
                                    fill_value=np.nan)
            else:
                logging.debug(
                    "No variables found for '" + str(VariableDefinition.LONG_NAME['_wet_binary_mask']) + "'")
                raise (VariableNameError("SymphonieReader",
                                         "No variables found for '" + str(
                                             VariableDefinition.LONG_NAME['wet_binary_mask']) + "'",
                                         1000))
        except Exception as ex:
            logging.debug("Error '" + str(ex) + "'")
            raise (VariableNameError("SymphonieReader", "An error occured : '" + str(ex) + "'", 1000))

    def read_variable_mesh_size(self, xmin, xmax, ymin, ymax):
        try:
            if "sqrt_dxdy" in self.grid.variables:
                return np.ma.filled(self.grid.variables["sqrt_dxdy"][ymin:ymax, xmin:xmax], fill_value=np.nan)
            else:
                logging.debug("No variables found for '" + str(VariableDefinition.LONG_NAME['mesh_size']) + "'")
                raise (VariableNameError("SymphonieReader",
                                         "No variables found for '" + str(
                                             VariableDefinition.LONG_NAME['mesh_size']) + "'",
                                         1000))
        except Exception as ex:
            logging.debug("Error '" + str(ex) + "'")
            raise (VariableNameError("SymphonieReader", "An error occured : '" + str(ex) + "'", 1000))

    def read_variable_x_mesh_size(self, xmin, xmax, ymin, ymax):
        try:
            if "dx_t" in self.grid.variables:
                return np.ma.filled(self.grid.variables["dx_t"][ymin:ymax, xmin:xmax], fill_value=np.nan)
            else:
                logging.debug("No variables found for '" + str(VariableDefinition.LONG_NAME['x_mesh_size']) + "'")
                raise (VariableNameError("SymphonieReader",
                                         "No variables found for '" + str(
                                             VariableDefinition.LONG_NAME['x_mesh_size']) + "'",
                                         1000))
        except Exception as ex:
            logging.debug("Error '" + str(ex) + "'")
            raise (VariableNameError("SymphonieReader", "An error occured : '" + str(ex) + "'", 1000))

    def read_variable_y_mesh_size(self, xmin, xmax, ymin, ymax):
        try:
            if "dy_t" in self.grid.variables:
                return np.ma.filled(self.grid.variables["dy_t"][ymin:ymax, xmin:xmax], fill_value=np.nan)
            else:
                logging.debug("No variables found for '" + str(VariableDefinition.LONG_NAME['y_mesh_size']) + "'")
                raise (VariableNameError("SymphonieReader",
                                         "No variables found for '" + str(
                                             VariableDefinition.LONG_NAME['y_mesh_size']) + "'",
                                         1000))
        except Exception as ex:
            logging.debug("Error '" + str(ex) + "'")
            raise (VariableNameError("SymphonieReader", "An error occured : '" + str(ex) + "'", 1000))

    #################
    # HYDRO
    # Sea Surface
    #################
    def read_variable_sea_surface_height_above_mean_sea_level_at_time(self, index_t, xmin, xmax, ymin, ymax):
        try:
            self.open_file(index_t)

            if "ssh_w" in self.ncfile.variables:
                data = np.ma.filled(self.ncfile.variables["ssh_w"][0, ymin:ymax, xmin:xmax], fill_value=np.nan)
            elif "ssh" in self.ncfile.variables:
                data = np.ma.filled(self.ncfile.variables["ssh"][0, ymin:ymax, xmin:xmax], fill_value=np.nan)
            elif "ssh_inst" in self.ncfile.variables:
                data = np.ma.filled(self.ncfile.variables["ssh_inst"][0, ymin:ymax, xmin:xmax], fill_value=np.nan)
            else:
                logging.debug("No variables found for '" + str(
                    VariableDefinition.LONG_NAME['sea_surface_height_above_mean_sea_level']) + "'")
                raise (VariableNameError("SymphonieReader",
                                         "No variables found for '" + str(
                                             VariableDefinition.LONG_NAME[
                                                 'sea_surface_height_above_mean_sea_level']) + "'",
                                         1000))

            if "wetmask_t" in self.ncfile.variables: # We apply the wetmask
                data[self.ncfile.variables["wetmask_t"][0, ymin:ymax, xmin:xmax] == 0] = np.nan

            return data

        except Exception as ex:
            logging.debug("Error '" + str(ex) + "'")
            raise (VariableNameError("SymphonieReader", "An error occured : '" + str(ex) + "'", 1000))

    def read_variable_sea_water_column_thickness_at_time(self, index_t, xmin, xmax, ymin, ymax):
        try:
            self.open_file(index_t)
            if "hssh" in self.ncfile.variables:
                data = np.ma.filled(self.ncfile.variables["hssh"][0, ymin:ymax, xmin:xmax], fill_value=np.nan)
            else:
                logging.debug("No variables found for '" + str(
                    VariableDefinition.LONG_NAME['sea_water_column_thickness']) + "'")
                raise (VariableNameError("SymphonieReader",
                                         "No variables found for '" + str(
                                             VariableDefinition.LONG_NAME['sea_water_column_thickness']) + "'",
                                         1000))

            if "wetmask_t" in self.ncfile.variables: # We apply the wetmask
                data[self.ncfile.variables["wetmask_t"][0, ymin:ymax, xmin:xmax] == 0] = np.nan

            return data

        except Exception as ex:
            logging.debug("Error '" + str(ex) + "'")
            raise (VariableNameError("SymphonieReader", "An error occured : '" + str(ex) + "'", 1000))

    def read_variable_sea_surface_temperature_at_time(self, index_t, xmin, xmax, ymin, ymax):
        try:
            self.open_file(index_t)
            index_z = self.get_z_size() - 1
            if "tem" in self.ncfile.variables:
                data = np.ma.filled(self.ncfile.variables["tem"][0, index_z, ymin:ymax, xmin:xmax],
                                    fill_value=np.nan)
            else:
                logging.debug("No variables found for '" + str(
                    VariableDefinition.LONG_NAME['sea_surface_temperature']) + "'")
                raise (VariableNameError("SymphonieReader",
                                         "No variables found for '" + str(
                                             VariableDefinition.LONG_NAME['sea_surface_temperature']) + "'",
                                         1000))

            if "wetmask_t" in self.ncfile.variables: # We apply the wetmask
                data[self.ncfile.variables["wetmask_t"][0, ymin:ymax, xmin:xmax] == 0] = np.nan

            return data

        except Exception as ex:
            logging.debug("Error '" + str(ex) + "'")
            raise (VariableNameError("SymphonieReader", "An error occured : '" + str(ex) + "'", 1000))



    def read_variable_sea_surface_salinity_at_time(self, index_t, xmin, xmax, ymin, ymax):
        try:
            self.open_file(index_t)
            index_z = self.get_z_size() - 1
            if "sal" in self.ncfile.variables:
                data =np.ma.filled(self.ncfile.variables["sal"][0, index_z, ymin:ymax, xmin:xmax],
                                    fill_value=np.nan)
            else:
                logging.debug("No variables found for '" + str(
                    VariableDefinition.LONG_NAME['sea_surface_salinity']) + "'")
                raise (VariableNameError("SymphonieReader",
                                         "No variables found for '" + str(
                                             VariableDefinition.LONG_NAME['sea_surface_salinity']) + "'",
                                         1000))

            if "wetmask_t" in self.ncfile.variables:  # We apply the wetmask
                data[self.ncfile.variables["wetmask_t"][0, ymin:ymax, xmin:xmax] == 0] = np.nan

            return data

        except Exception as ex:
            logging.debug("Error '" + str(ex) + "'")
            raise (VariableNameError("SymphonieReader", "An error occured : '" + str(ex) + "'", 1000))

    def read_variable_sea_water_velocity_at_sea_water_surface_at_time(self, index_t, xmin, xmax, ymin, ymax):
        try:
            self.open_file(index_t)
            index_z = self.get_z_size() - 1
            xmin_overlap, xmax_overlap, ymin_overlap, ymax_overlap, new_xmin, new_xmax, new_ymin, new_ymax = self.compute_overlap_indexes(
                xmin, xmax, ymin, ymax)

            mask_t = self.grid.variables["mask_t"][index_z, ymin_overlap:ymax_overlap, xmin_overlap:xmax_overlap];
            mask_u = self.grid.variables["mask_u"][index_z, ymin_overlap:ymax_overlap, xmin_overlap:xmax_overlap];
            mask_v = self.grid.variables["mask_v"][index_z, ymin_overlap:ymax_overlap, xmin_overlap:xmax_overlap];

            if self.gridrotcos_t is None and self.gridrotsin_t is None:
                self.compute_rot()

            rotcos = self.gridrotcos_t[ymin_overlap:ymax_overlap, xmin_overlap:xmax_overlap];
            rotsin = self.gridrotsin_t[ymin_overlap:ymax_overlap, xmin_overlap:xmax_overlap];

            if "vel_u" in self.ncfile.variables and "vel_v" in self.ncfile.variables:
                data_u = np.ma.filled(
                    self.ncfile.variables["vel_u"][0, index_z, ymin_overlap:ymax_overlap, xmin_overlap:xmax_overlap],
                    fill_value=np.nan)
                data_v = np.ma.filled(
                    self.ncfile.variables["vel_v"][0, index_z, ymin_overlap:ymax_overlap, xmin_overlap:xmax_overlap],
                    fill_value=np.nan)
            elif "u" in self.ncfile.variables and "v" in self.ncfile.variables:
                data_u = np.ma.filled(
                    self.ncfile.variables["u"][0, index_z, ymin_overlap:ymax_overlap, xmin_overlap:xmax_overlap],
                    fill_value=np.nan)
                data_v = np.ma.filled(
                    self.ncfile.variables["v"][0, index_z, ymin_overlap:ymax_overlap, xmin_overlap:xmax_overlap],
                    fill_value=np.nan)
            else:
                logging.debug("No variables found for 'Sea Water Velocity At Sea Water Surface'")
                raise (VariableNameError("SymphonieReader",
                                         "No variables found for 'Sea Water Velocity At Sea Water Surface'",
                                         1000))

            u_rot, v_rot = self.compute_vector_rotation(data_u, data_v, rotcos, rotsin, mask_t, mask_u, mask_v)

            if "wetmask_t" in self.ncfile.variables: # We apply the wetmask
                u_rot[self.ncfile.variables["wetmask_t"][0, ymin_overlap:ymax_overlap,
                      xmin_overlap:xmax_overlap] == 0] = np.nan
                v_rot[self.ncfile.variables["wetmask_t"][0, ymin_overlap:ymax_overlap,
                      xmin_overlap:xmax_overlap] == 0] = np.nan

            return [u_rot[new_ymin:new_ymax, new_xmin:new_xmax], v_rot[new_ymin:new_ymax, new_xmin:new_xmax]]

        except Exception as ex:
            logging.debug("Error '" + str(ex) + "'")
            raise (VariableNameError("SymphonieReader", "An error occured : '" + str(ex) + "'", 1000))

    #################
    # HYDRO
    # Ground level
    #################

    def read_variable_sea_water_temperature_at_ground_level_at_time(self, index_t, xmin, xmax, ymin, ymax):
        try:
            self.open_file(index_t)
            index_z = 0
            if "tem" in self.ncfile.variables:
                data = np.ma.filled(self.ncfile.variables["tem"][0, index_z, ymin:ymax, xmin:xmax],
                                    fill_value=np.nan)
            else:
                logging.debug("No variables found for '" + str(
                    VariableDefinition.LONG_NAME['sea_water_temperature_at_ground_level']) + "'")
                raise (VariableNameError("SymphonieReader",
                                         "No variables found for '" + str(
                                             VariableDefinition.LONG_NAME['sea_water_temperature_at_ground_level']) + "'",
                                         1000))

            if "wetmask_t" in self.ncfile.variables: # We apply the wetmask
                data[self.ncfile.variables["wetmask_t"][0, ymin:ymax, xmin:xmax] == 0] = np.nan

            return data

        except Exception as ex:
            logging.debug("Error '" + str(ex) + "'")
            raise (VariableNameError("SymphonieReader", "An error occured : '" + str(ex) + "'", 1000))

    def read_variable_sea_water_salinity_at_ground_level_at_time(self, index_t, xmin, xmax, ymin, ymax):
        try:
            self.open_file(index_t)
            index_z = 0
            if "sal" in self.ncfile.variables:
                data = np.ma.filled(self.ncfile.variables["sal"][0, index_z, ymin:ymax, xmin:xmax],
                                    fill_value=np.nan)
            else:
                logging.debug("No variables found for '" + str(
                    VariableDefinition.LONG_NAME['sea_water_salinity_at_ground_level']) + "'")
                raise (VariableNameError("SymphonieReader",
                                         "No variables found for '" + str(
                                             VariableDefinition.LONG_NAME['sea_water_salinity_at_ground_level']) + "'",
                                         1000))

            if "wetmask_t" in self.ncfile.variables: # We apply the wetmask
                data[self.ncfile.variables["wetmask_t"][0, ymin:ymax, xmin:xmax] == 0] = np.nan

            return data

        except Exception as ex:
            logging.debug("Error '" + str(ex) + "'")
            raise (VariableNameError("SymphonieReader", "An error occured : '" + str(ex) + "'", 1000))


    def read_variable_sea_water_velocity_at_ground_level_at_time(self, index_t, xmin, xmax, ymin, ymax):
        try:
            self.open_file(index_t)
            index_z = 0
            xmin_overlap, xmax_overlap, ymin_overlap, ymax_overlap, new_xmin, new_xmax, new_ymin, new_ymax = self.compute_overlap_indexes(
                xmin, xmax, ymin, ymax)

            mask_t = self.grid.variables["mask_t"][index_z, ymin_overlap:ymax_overlap, xmin_overlap:xmax_overlap];
            mask_u = self.grid.variables["mask_u"][index_z, ymin_overlap:ymax_overlap, xmin_overlap:xmax_overlap];
            mask_v = self.grid.variables["mask_v"][index_z, ymin_overlap:ymax_overlap, xmin_overlap:xmax_overlap];

            if self.gridrotcos_t is None and self.gridrotsin_t is None:
                self.compute_rot()

            rotcos = self.gridrotcos_t[ymin_overlap:ymax_overlap, xmin_overlap:xmax_overlap];
            rotsin = self.gridrotsin_t[ymin_overlap:ymax_overlap, xmin_overlap:xmax_overlap];

            if "vel_u" in self.ncfile.variables and "vel_v" in self.ncfile.variables:
                data_u = np.ma.filled(
                    self.ncfile.variables["vel_u"][0, index_z, ymin_overlap:ymax_overlap,
                    xmin_overlap:xmax_overlap],
                    fill_value=np.nan)
                data_v = np.ma.filled(
                    self.ncfile.variables["vel_v"][0, index_z, ymin_overlap:ymax_overlap,
                    xmin_overlap:xmax_overlap],
                    fill_value=np.nan)
            elif "u" in self.ncfile.variables and "v" in self.ncfile.variables:
                data_u = np.ma.filled(
                    self.ncfile.variables["u"][0, index_z, ymin_overlap:ymax_overlap,
                    xmin_overlap:xmax_overlap],
                    fill_value=np.nan)
                data_v = np.ma.filled(
                    self.ncfile.variables["v"][0, index_z, ymin_overlap:ymax_overlap,
                    xmin_overlap:xmax_overlap],
                    fill_value=np.nan)
            else:
                logging.debug("No variables found for 'Sea Water Velocity at Ground level'")
                raise (VariableNameError("SymphonieReader",
                                         "No variables found for 'Sea Water Velocity at Ground level'",
                                         1000))

            u_rot, v_rot = self.compute_vector_rotation(data_u, data_v, rotcos, rotsin, mask_t, mask_u, mask_v)

            if "wetmask_t" in self.ncfile.variables: # We apply the wetmask
                u_rot[self.ncfile.variables["wetmask_t"][0, ymin_overlap:ymax_overlap,
                      xmin_overlap:xmax_overlap] == 0] = np.nan
                v_rot[self.ncfile.variables["wetmask_t"][0, ymin_overlap:ymax_overlap,
                      xmin_overlap:xmax_overlap] == 0] = np.nan

            return [u_rot[new_ymin:new_ymax, new_xmin:new_xmax], v_rot[new_ymin:new_ymax, new_xmin:new_xmax]]

        except Exception as ex:
            logging.debug("Error '" + str(ex) + "'")
            raise (VariableNameError("SymphonieReader", "An error occured : '" + str(ex) + "'", 1000))

    #################
    # HYDRO
    # 2D
    #################

    def read_variable_bathymetry(self, xmin, xmax, ymin, ymax):
        try:
            if "hm_w" in self.grid.variables:
                return np.ma.filled(self.grid.variables["hm_w"][ymin:ymax, xmin:xmax], fill_value=np.nan)
            elif "h_w" in self.grid.variables:
                return np.ma.filled(self.grid.variables["h_w"][ymin:ymax, xmin:xmax], fill_value=np.nan)
            else:
                logging.debug("No variables found for '" + str(VariableDefinition.LONG_NAME['bathymetry']) + "'")
                raise (VariableNameError("SymphonieReader",
                                         "No variables found for '" + str(
                                             VariableDefinition.LONG_NAME['bathymetry']) + "'",
                                         1000))

        except Exception as ex:
            logging.debug("Error '" + str(ex) + "'")
            raise (VariableNameError("SymphonieReader", "An error occured : '" + str(ex) + "'", 1000))

    def read_variable_barotropic_sea_water_velocity_at_time(self, index_t, xmin, xmax, ymin, ymax):
        try:
            self.open_file(index_t)
            index_z = self.get_z_size() - 1
            xmin_overlap, xmax_overlap, ymin_overlap, ymax_overlap, new_xmin, new_xmax, new_ymin, new_ymax = self.compute_overlap_indexes(
                xmin, xmax, ymin, ymax)

            mask_t = self.grid.variables["mask_t"][index_z, ymin_overlap:ymax_overlap, xmin_overlap:xmax_overlap];
            mask_u = self.grid.variables["mask_u"][index_z, ymin_overlap:ymax_overlap, xmin_overlap:xmax_overlap];
            mask_v = self.grid.variables["mask_v"][index_z, ymin_overlap:ymax_overlap, xmin_overlap:xmax_overlap];

            if self.gridrotcos_t is None and self.gridrotsin_t is None:
                self.compute_rot()

            rotcos = self.gridrotcos_t[ymin_overlap:ymax_overlap, xmin_overlap:xmax_overlap];
            rotsin = self.gridrotsin_t[ymin_overlap:ymax_overlap, xmin_overlap:xmax_overlap];

            if "velbar_u" in self.ncfile.variables and "velbar_v" in self.ncfile.variables:
                data_u = np.ma.filled(
                    self.ncfile.variables["velbar_u"][0, ymin_overlap:ymax_overlap, xmin_overlap:xmax_overlap],
                    fill_value=np.nan)
                data_v = np.ma.filled(
                    self.ncfile.variables["velbar_v"][0, ymin_overlap:ymax_overlap, xmin_overlap:xmax_overlap],
                    fill_value=np.nan)
            else:
                logging.debug("No variables found for 'Barotropic Sea Water Velocity'")
                raise (VariableNameError("SymphonieReader",
                                         "No variables found for 'Barotropic Sea Water Velocity'",
                                         1000))

            u_rot, v_rot = self.compute_vector_rotation(data_u, data_v, rotcos, rotsin, mask_t, mask_u, mask_v)

            if "wetmask_t" in self.ncfile.variables:  # We apply the wetmask
                u_rot[self.ncfile.variables["wetmask_t"][0, ymin_overlap:ymax_overlap,
                      xmin_overlap:xmax_overlap] == 0] = np.nan
                v_rot[self.ncfile.variables["wetmask_t"][0, ymin_overlap:ymax_overlap,
                      xmin_overlap:xmax_overlap] == 0] = np.nan

            return [u_rot[new_ymin:new_ymax, new_xmin:new_xmax], v_rot[new_ymin:new_ymax, new_xmin:new_xmax]]

        except Exception as ex:
            logging.debug("Error '" + str(ex) + "'")
            raise (VariableNameError("SymphonieReader", "An error occured : '" + str(ex) + "'", 1000))

    #################
    # HYDRO
    # 3D
    #################

    def read_variable_depth_at_depth(self, index_z, xmin, xmax, ymin, ymax):
        try:
            if "depth_t" in self.grid.variables:
                data = self.grid.variables["depth_t"][index_z, ymin:ymax, xmin:xmax]
                data[::] *= -1.0  # inverse la profondeur
                return np.ma.filled(data, fill_value=np.nan)
            else:
                logging.debug("No variables found for '" + str(VariableDefinition.LONG_NAME['depth_sigma']) + "'")
                raise (VariableNameError("SymphonieReader",
                                         "No variables found for '" + str(
                                             VariableDefinition.LONG_NAME['depth_sigma']) + "'",
                                         1000))
        except Exception as ex:
            logging.debug("Error '" + str(ex) + "'")
            raise (VariableNameError("SymphonieReader", "An error occured : '" + str(ex) + "'", 1000))

    def read_variable_sea_water_temperature_at_time_and_depth(self, index_t, index_z, xmin, xmax, ymin, ymax):
        try:
            self.open_file(index_t)
            if "tem" in self.ncfile.variables:
                data = np.ma.filled(self.ncfile.variables["tem"][0, index_z, ymin:ymax, xmin:xmax],
                                    fill_value=np.nan)
            else:
                logging.debug("No variables found for '" + str(
                    VariableDefinition.LONG_NAME['sea_water_temperature']) + "'")
                raise (VariableNameError("SymphonieReader",
                                         "No variables found for '" + str(
                                             VariableDefinition.LONG_NAME['sea_water_temperature']) + "'",
                                         1000))

            if "wetmask_t" in self.ncfile.variables:  # We apply the wetmask
                data[self.ncfile.variables["wetmask_t"][0, ymin:ymax, xmin:xmax] == 0] = np.nan

            return data

        except Exception as ex:
            logging.debug("Error '" + str(ex) + "'")
            raise (VariableNameError("SymphonieReader", "An error occured : '" + str(ex) + "'", 1000))

    def read_variable_sea_water_salinity_at_time_and_depth(self, index_t, index_z, xmin, xmax, ymin, ymax):
        try:
            self.open_file(index_t)
            if "sal" in self.ncfile.variables:
                data = np.ma.filled(self.ncfile.variables["sal"][0, index_z, ymin:ymax, xmin:xmax],
                                    fill_value=np.nan)
            else:
                logging.debug("No variables found for '" + str(
                    VariableDefinition.LONG_NAME['sea_water_salinity']) + "'")
                raise (VariableNameError("SymphonieReader",
                                         "No variables found for '" + str(
                                             VariableDefinition.LONG_NAME['sea_water_salinity']) + "'",
                                         1000))

            if "wetmask_t" in self.ncfile.variables: # We apply the wetmask
                data[self.ncfile.variables["wetmask_t"][0, ymin:ymax, xmin:xmax] == 0] = np.nan

            return data

        except Exception as ex:
            logging.debug("Error '" + str(ex) + "'")
            raise (VariableNameError("SymphonieReader", "An error occured : '" + str(ex) + "'", 1000))

    def read_variable_baroclinic_sea_water_velocity_at_time_and_depth(self, index_t, index_z, xmin, xmax, ymin, ymax):
        try:
            self.open_file(index_t)
            xmin_overlap, xmax_overlap, ymin_overlap, ymax_overlap, new_xmin, new_xmax, new_ymin, new_ymax = self.compute_overlap_indexes(
                xmin, xmax, ymin, ymax)

            mask_t = self.grid.variables["mask_t"][index_z, ymin_overlap:ymax_overlap, xmin_overlap:xmax_overlap];
            mask_u = self.grid.variables["mask_u"][index_z, ymin_overlap:ymax_overlap, xmin_overlap:xmax_overlap];
            mask_v = self.grid.variables["mask_v"][index_z, ymin_overlap:ymax_overlap, xmin_overlap:xmax_overlap];

            if self.gridrotcos_t is None and self.gridrotsin_t is None:
                self.compute_rot()

            rotcos = self.gridrotcos_t[ymin_overlap:ymax_overlap, xmin_overlap:xmax_overlap];
            rotsin = self.gridrotsin_t[ymin_overlap:ymax_overlap, xmin_overlap:xmax_overlap];

            if "vel_u" in self.ncfile.variables and  "vel_v" in self.ncfile.variables:
                data_u = np.ma.filled(
                    self.ncfile.variables["vel_u"][0, index_z, ymin_overlap:ymax_overlap,
                    xmin_overlap:xmax_overlap],
                    fill_value=np.nan)
                data_v = np.ma.filled(
                    self.ncfile.variables["vel_v"][0, index_z, ymin_overlap:ymax_overlap,
                    xmin_overlap:xmax_overlap],
                    fill_value=np.nan)
            elif "u" in self.ncfile.variables and "v" in self.ncfile.variables:
                data_u = np.ma.filled(
                    self.ncfile.variables["u"][0, index_z, ymin_overlap:ymax_overlap,
                    xmin_overlap:xmax_overlap],
                    fill_value=np.nan)
                data_v = np.ma.filled(
                    self.ncfile.variables["v"][0, index_z, ymin_overlap:ymax_overlap,
                    xmin_overlap:xmax_overlap],
                    fill_value=np.nan)
            else:
                logging.debug("No variables found for 'Baroclinic Sea Water Velocity'")
                raise (VariableNameError("SymphonieReader",
                                         "No variables found for 'Baroclinic Sea Water Velocity'",
                                         1000))

            u_rot, v_rot = self.compute_vector_rotation(data_u, data_v, rotcos, rotsin, mask_t, mask_u, mask_v)

            if "wetmask_t" in self.ncfile.variables: # We apply the wetmask
                u_rot[self.ncfile.variables["wetmask_t"][0, ymin_overlap:ymax_overlap,
                      xmin_overlap:xmax_overlap] == 0] = np.nan
                v_rot[self.ncfile.variables["wetmask_t"][0, ymin_overlap:ymax_overlap,
                      xmin_overlap:xmax_overlap] == 0] = np.nan

            return [u_rot[new_ymin:new_ymax, new_xmin:new_xmax], v_rot[new_ymin:new_ymax, new_xmin:new_xmax]]

        except Exception as ex:
            logging.debug("Error '" + str(ex) + "'")
            raise (VariableNameError("SymphonieReader", "An error occured : '" + str(ex) + "'", 1000))

    #################
    # WAVES
    # Sea Surface
    #################
    def read_variable_sea_surface_wave_significant_height_at_time(self, index_t, xmin, xmax, ymin, ymax):
        try:
            self.open_file(index_t)
            if "hs_wave_t" in self.ncfile.variables:
                data = np.ma.filled(self.ncfile.variables["hs_wave_t"][0, ymin:ymax, xmin:xmax],
                                    fill_value=np.nan)
            else:
                logging.debug(
                    "No variables found for '" + str(
                        VariableDefinition.LONG_NAME['sea_surface_wave_significant_height']) + "'")
                raise (VariableNameError("SymphonieReader",
                                         "No variables found for '" + str(
                                             VariableDefinition.LONG_NAME['sea_surface_wave_significant_height']) + "'",
                                         1000))

            if "wetmask_t" in self.ncfile.variables: # We apply the wetmask
                data[self.ncfile.variables["wetmask_t"][0, ymin:ymax, xmin:xmax] == 0] = np.nan

            return data

        except Exception as ex:
            logging.debug("Error '" + str(ex) + "'")
            raise (VariableNameError("SymphonieReader", "An error occured : '" + str(ex) + "'", 1000))

    def read_variable_sea_surface_wave_mean_period_at_time(self, index_t, xmin, xmax, ymin, ymax):
        try:
            self.open_file(index_t)
            if "t_wave_t" in self.ncfile.variables:
                data = np.ma.filled(self.ncfile.variables["t_wave_t"][0, ymin:ymax, xmin:xmax], fill_value=np.nan)
            else:
                logging.debug(
                    "No variables found for '" + str(
                        VariableDefinition.LONG_NAME['sea_surface_wave_mean_period']) + "'")
                raise (VariableNameError("SymphonieReader",
                                         "No variables found for '" + str(
                                             VariableDefinition.LONG_NAME['sea_surface_wave_mean_period']) + "'",
                                         1000))

            if "wetmask_t" in self.ncfile.variables: # We apply the wetmask
                data[self.ncfile.variables["wetmask_t"][0, ymin:ymax, xmin:xmax] == 0] = np.nan

            return data

        except Exception as ex:
            logging.debug("Error '" + str(ex) + "'")
            raise (VariableNameError("SymphonieReader", "An error occured : '" + str(ex) + "'", 1000))

    def read_variable_sea_surface_wave_to_direction_at_time(self, index_t, xmin, xmax, ymin, ymax):
        try:
            self.open_file(index_t)
            if "dir_wave_t" in self.ncfile.variables:
                data = np.ma.filled(self.ncfile.variables["dir_wate_t"][0, ymin:ymax, xmin:xmax],
                                    fill_value=np.nan)
            else:
                logging.debug(
                    "No variables found for '" + str(
                        VariableDefinition.LONG_NAME['sea_surface_wave_to_direction']) + "'")
                raise (VariableNameError("SymphonieReader",
                                         "No variables found for '" + str(
                                             VariableDefinition.LONG_NAME['sea_surface_wave_to_direction']) + "'",
                                         1000))

            if "wetmask_t" in self.ncfile.variables: # We apply the wetmask
                data[self.ncfile.variables["wetmask_t"][0, ymin:ymax, xmin:xmax] == 0] = np.nan

            return data

        except Exception as ex:
            logging.debug("Error '" + str(ex) + "'")
            raise (VariableNameError("SymphonieReader", "An error occured : '" + str(ex) + "'", 1000))

    def read_variable_sea_surface_wave_stokes_drift_velocity_at_time(self, index_t, xmin, xmax, ymin, ymax):
        try:
            self.open_file(index_t)
            index_z = self.get_z_size() - 1
            xmin_overlap, xmax_overlap, ymin_overlap, ymax_overlap, new_xmin, new_xmax, new_ymin, new_ymax = self.compute_overlap_indexes(
                xmin, xmax, ymin, ymax)

            mask_t = self.grid.variables["mask_t"][index_z, ymin_overlap:ymax_overlap, xmin_overlap:xmax_overlap];
            mask_u = self.grid.variables["mask_u"][index_z, ymin_overlap:ymax_overlap, xmin_overlap:xmax_overlap];
            mask_v = self.grid.variables["mask_v"][index_z, ymin_overlap:ymax_overlap, xmin_overlap:xmax_overlap];

            if self.gridrotcos_t is None and self.gridrotsin_t is None:
                self.compute_rot()

            rotcos = self.gridrotcos_t[ymin_overlap:ymax_overlap, xmin_overlap:xmax_overlap];
            rotsin = self.gridrotsin_t[ymin_overlap:ymax_overlap, xmin_overlap:xmax_overlap];

            if "velbarstokes_u" in self.ncfile.variables and "velbarstokes_v" in self.ncfile.variables:
                data_u = np.ma.filled(
                    self.ncfile.variables["velbarstokes_u"][0, ymin_overlap:ymax_overlap,
                    xmin_overlap:xmax_overlap],
                    fill_value=np.nan)
                data_v = np.ma.filled(
                    self.ncfile.variables["velbarstokes_v"][0, ymin_overlap:ymax_overlap,
                    xmin_overlap:xmax_overlap],
                    fill_value=np.nan)
            else:
                logging.debug("No variables found for 'Surface Stokes Drift Velocity'")
                raise (VariableNameError("SymphonieReader",
                                         "No variables found for 'Surface Stokes Drift Velocity'",
                                         1000))

            u_rot, v_rot = self.compute_vector_rotation(data_u, data_v, rotcos, rotsin, mask_t, mask_u, mask_v)

            if "wetmask_t" in self.ncfile.variables: # We apply the wetmask
                u_rot[self.ncfile.variables["wetmask_t"][0, ymin_overlap:ymax_overlap,
                      xmin_overlap:xmax_overlap] == 0] = np.nan
                v_rot[self.ncfile.variables["wetmask_t"][0, ymin_overlap:ymax_overlap,
                      xmin_overlap:xmax_overlap] == 0] = np.nan

            return [u_rot[new_ymin:new_ymax, new_xmin:new_xmax], v_rot[new_ymin:new_ymax, new_xmin:new_xmax]]

        except Exception as ex:
            logging.debug("Error '" + str(ex) + "'")
            raise (VariableNameError("SymphonieReader", "An error occured : '" + str(ex) + "'", 1000))

    #################
    # WAVES
    # Momentum flux
    #################

    def read_variable_atmosphere_momentum_flux_to_waves_at_time(self, index_t, xmin, xmax, ymin, ymax):
        try:
            self.open_file(index_t)
            index_z = self.get_z_size() - 1
            xmin_overlap, xmax_overlap, ymin_overlap, ymax_overlap, new_xmin, new_xmax, new_ymin, new_ymax = self.compute_overlap_indexes(
                xmin, xmax, ymin, ymax)

            mask_t = self.grid.variables["mask_t"][index_z, ymin_overlap:ymax_overlap, xmin_overlap:xmax_overlap];
            mask_u = self.grid.variables["mask_u"][index_z, ymin_overlap:ymax_overlap, xmin_overlap:xmax_overlap];
            mask_v = self.grid.variables["mask_v"][index_z, ymin_overlap:ymax_overlap, xmin_overlap:xmax_overlap];

            if self.gridrotcos_t is None and self.gridrotsin_t is None:
                self.compute_rot()

            rotcos = self.gridrotcos_t[ymin_overlap:ymax_overlap, xmin_overlap:xmax_overlap];
            rotsin = self.gridrotsin_t[ymin_overlap:ymax_overlap, xmin_overlap:xmax_overlap];

            if "tawx" in self.ncfile.variables and  "tawy" in self.ncfile.variables: # We apply the wetmask
                data_u = np.ma.filled(
                    self.ncfile.variables["tawx"][0, ymin_overlap:ymax_overlap, xmin_overlap:xmax_overlap],
                    fill_value=np.nan)
                data_v = np.ma.filled(
                    self.ncfile.variables["tawy"][0, ymin_overlap:ymax_overlap, xmin_overlap:xmax_overlap],
                    fill_value=np.nan)
            else:
                logging.debug("No variables found for 'Atmosphere Momentum Flux to Waves'")
                raise (VariableNameError("SymphonieReader",
                                         "No variables found for 'Atmosphere Momentum Flux to Waves'",
                                         1000))

            u_rot, v_rot = self.compute_vector_rotation(data_u, data_v, rotcos, rotsin, mask_t, mask_u, mask_v)

            if "wetmask_t" in self.ncfile.variables: # We apply the wetmask
                u_rot[self.ncfile.variables["wetmask_t"][0, ymin_overlap:ymax_overlap,
                      xmin_overlap:xmax_overlap] == 0] = np.nan
                v_rot[self.ncfile.variables["wetmask_t"][0, ymin_overlap:ymax_overlap,
                      xmin_overlap:xmax_overlap] == 0] = np.nan

            return [u_rot[new_ymin:new_ymax, new_xmin:new_xmax], v_rot[new_ymin:new_ymax, new_xmin:new_xmax]]

        except Exception as ex:
            logging.debug("Error '" + str(ex) + "'")
            raise (VariableNameError("SymphonieReader", "An error occured : '" + str(ex) + "'", 1000))

    def read_variable_waves_momentum_flux_to_ocean_at_time(self, index_t, xmin, xmax, ymin, ymax):
        try:
            self.open_file(index_t)
            index_z = self.get_z_size() - 1
            xmin_overlap, xmax_overlap, ymin_overlap, ymax_overlap, new_xmin, new_xmax, new_ymin, new_ymax = self.compute_overlap_indexes(
                xmin, xmax, ymin, ymax)

            mask_t = self.grid.variables["mask_t"][index_z, ymin_overlap:ymax_overlap, xmin_overlap:xmax_overlap];
            mask_u = self.grid.variables["mask_u"][index_z, ymin_overlap:ymax_overlap, xmin_overlap:xmax_overlap];
            mask_v = self.grid.variables["mask_v"][index_z, ymin_overlap:ymax_overlap, xmin_overlap:xmax_overlap];

            if self.gridrotcos_t is None and self.gridrotsin_t is None:
                self.compute_rot()

            rotcos = self.gridrotcos_t[ymin_overlap:ymax_overlap, xmin_overlap:xmax_overlap];
            rotsin = self.gridrotsin_t[ymin_overlap:ymax_overlap, xmin_overlap:xmax_overlap];

            if "twox" in self.ncfile.variables and "twoy" in self.ncfile.variables:
                data_u = np.ma.filled(
                    self.ncfile.variables["twox"][0, ymin_overlap:ymax_overlap, xmin_overlap:xmax_overlap],
                    fill_value=np.nan)
                data_v = np.ma.filled(
                    self.ncfile.variables["twoy"][0, ymin_overlap:ymax_overlap, xmin_overlap:xmax_overlap],
                    fill_value=np.nan)
            else:
                logging.debug("No variables found for 'Waves Momentum Flux To Ocean'")
                raise (VariableNameError("SymphonieReader",
                                         "No variables found for 'Waves Momentum Flux To Ocean'",
                                         1000))

            u_rot, v_rot = self.compute_vector_rotation(data_u, data_v, rotcos, rotsin, mask_t, mask_u, mask_v)

            if "wetmask_t" in self.ncfile.variables: # We apply the wetmask
                u_rot[self.ncfile.variables["wetmask_t"][0, ymin_overlap:ymax_overlap,
                      xmin_overlap:xmax_overlap] == 0] = np.nan
                v_rot[self.ncfile.variables["wetmask_t"][0, ymin_overlap:ymax_overlap,
                      xmin_overlap:xmax_overlap] == 0] = np.nan

            return [u_rot[new_ymin:new_ymax, new_xmin:new_xmax], v_rot[new_ymin:new_ymax, new_xmin:new_xmax]]

        except Exception as ex:
            logging.debug("Error '" + str(ex) + "'")
            raise (VariableNameError("SymphonieReader", "An error occured : '" + str(ex) + "'", 1000))

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

            mask_t = self.grid.variables["mask_t"][index_z, ymin_overlap:ymax_overlap, xmin_overlap:xmax_overlap];
            mask_u = self.grid.variables["mask_u"][index_z, ymin_overlap:ymax_overlap, xmin_overlap:xmax_overlap];
            mask_v = self.grid.variables["mask_v"][index_z, ymin_overlap:ymax_overlap, xmin_overlap:xmax_overlap];

            if self.gridrotcos_t is None and self.gridrotsin_t is None:
                self.compute_rot()

            rotcos = self.gridrotcos_t[ymin_overlap:ymax_overlap, xmin_overlap:xmax_overlap];
            rotsin = self.gridrotsin_t[ymin_overlap:ymax_overlap, xmin_overlap:xmax_overlap];

            if "wstress_u" in self.ncfile.variables and "wstress_v" in self.ncfile.variables:
                data_u = np.ma.filled(
                    self.ncfile.variables["wstress_u"][0, ymin_overlap:ymax_overlap, xmin_overlap:xmax_overlap],
                    fill_value=np.nan)
                data_v = np.ma.filled(
                    self.ncfile.variables["wstress_v"][0, ymin_overlap:ymax_overlap, xmin_overlap:xmax_overlap],
                    fill_value=np.nan)
            else:
                logging.debug("No variables found for 'Wind Stress'")
                raise (VariableNameError("SymphonieReader",
                                         "No variables found for 'Wind Stress'",
                                         1000))

            u_rot, v_rot = self.compute_vector_rotation(data_u, data_v, rotcos, rotsin, mask_t, mask_u, mask_v)

            if "wetmask_t" in self.ncfile.variables: # We apply the wetmask
                u_rot[self.ncfile.variables["wetmask_t"][0, ymin_overlap:ymax_overlap,
                      xmin_overlap:xmax_overlap] == 0] = np.nan
                v_rot[self.ncfile.variables["wetmask_t"][0, ymin_overlap:ymax_overlap,
                      xmin_overlap:xmax_overlap] == 0] = np.nan

            return [u_rot[new_ymin:new_ymax, new_xmin:new_xmax], v_rot[new_ymin:new_ymax, new_xmin:new_xmax]]

        except Exception as ex:
            logging.debug("Error '" + str(ex) + "'")
            raise (VariableNameError("SymphonieReader", "An error occured : '" + str(ex) + "'", 1000))

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

            mask_t = self.grid.variables["mask_t"][index_z, ymin_overlap:ymax_overlap, xmin_overlap:xmax_overlap];
            mask_u = self.grid.variables["mask_u"][index_z, ymin_overlap:ymax_overlap, xmin_overlap:xmax_overlap];
            mask_v = self.grid.variables["mask_v"][index_z, ymin_overlap:ymax_overlap, xmin_overlap:xmax_overlap];

            if self.gridrotcos_t is None and self.gridrotsin_t is None:
                self.compute_rot()

            rotcos = self.gridrotcos_t[ymin_overlap:ymax_overlap, xmin_overlap:xmax_overlap];
            rotsin = self.gridrotsin_t[ymin_overlap:ymax_overlap, xmin_overlap:xmax_overlap];

            if "uwind_t" in self.ncfile.variables and "vwind_t" in self.ncfile.variables:
                data_u = np.ma.filled(
                    self.ncfile.variables["uwind_t"][0, ymin_overlap:ymax_overlap, xmin_overlap:xmax_overlap],
                    fill_value=np.nan)
                data_v = np.ma.filled(
                    self.ncfile.variables["vwind_t"][0, ymin_overlap:ymax_overlap, xmin_overlap:xmax_overlap],
                    fill_value=np.nan)
            else:
                logging.debug("No variables found for \'Wind 10m\'")
                raise (VariableNameError("SymphonieReader",
                                         "No variables found for \'Wind 10m\'",
                                         1000))

            u_rot, v_rot = self.compute_vector_rotation(data_u, data_v, rotcos, rotsin, mask_t, mask_u, mask_v)

            if "wetmask_t" in self.ncfile.variables: # We apply the wetmask
                u_rot[self.ncfile.variables["wetmask_t"][0, ymin_overlap:ymax_overlap,
                      xmin_overlap:xmax_overlap] == 0] = np.nan
                v_rot[self.ncfile.variables["wetmask_t"][0, ymin_overlap:ymax_overlap,
                      xmin_overlap:xmax_overlap] == 0] = np.nan

            return [u_rot[new_ymin:new_ymax, new_xmin:new_xmax], v_rot[new_ymin:new_ymax, new_xmin:new_xmax]]

        except Exception as ex:
            logging.debug("Error '" + str(ex) + "'")
            raise (VariableNameError("SymphonieReader", "An error occured : '" + str(ex) + "'", 1000))

    # Others
    def read_variable_Ha(self, xmin, xmax, ymin, ymax):
        try:
            if "Ha" in self.ncfile.variables:
                data = np.ma.filled(self.ncfile.variables["Ha"][ymin:ymax, xmin:xmax],
                                    fill_value=np.nan)
            else:
                logging.debug(
                    "No variables found for 'Ha'")
                raise (VariableNameError("SymphonieReader",
                                         "No variables found for 'Ha'",
                                         1000))

            if "wetmask_t" in self.ncfile.variables: # We apply the wetmask
                data[self.ncfile.variables["wetmask_t"][0, ymin:ymax, xmin:xmax] == 0] = np.nan

            return data

        except Exception as ex:
            logging.debug("Error '" + str(ex) + "'")
            raise (VariableNameError("SymphonieReader", "An error occured : '" + str(ex) + "'", 1000))
