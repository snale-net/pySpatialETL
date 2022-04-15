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

from datetime import datetime
from datetime import timedelta

import cfgrib
import numpy as np

from spatialetl.coverage.TimeCoverage import TimeCoverage
from spatialetl.coverage.io.CoverageReader import CoverageReader
from spatialetl.exception.VariableNameError import VariableNameError
from spatialetl.utils.VariableDefinition import VariableDefinition
from spatialetl.utils.logger import logging


class ERA5Reader (CoverageReader):

    def __init__(self, myFile):
        CoverageReader.__init__(self,myFile);

        try:
            self.ds_instant = cfgrib.open_file(self.filename, filter_by_keys={'stepType': 'instant'})
        except Exception as ex:
            logging.debug("Error '" + str(ex) + "'")
            try:
                self.ds_instant = cfgrib.open_file(self.filename, filter_by_keys={'stepType': 'accum'})
            except Exception as ex:
                logging.debug("Error '" + str(ex) + "'")
                raise VariableNameError("ECMWFReader", "An error occured : '" + str(ex) + "'", 1000)

        self.ds_accum = cfgrib.open_file(self.filename, filter_by_keys={'stepType': 'accum'})

        logging.debug(sorted(self.ds_instant.variables))
        logging.debug(sorted(self.ds_accum.variables))

        lon = self.ds_instant.variables['longitude'].data
        self.new_lon = np.mod(lon + 180.0, 360.0) - 180.0
        lat = self.ds_instant.variables['latitude'].data
        xx, yy = np.meshgrid(self.new_lon, lat)
        self.new_indexes = np.argsort(xx)

        # if os.path.isfile(self.filename):
        #     self.ncfile = Dataset(self.filename, 'r')
        # elif os.path.isdir(self.filename):
        #     self.ncfile = MFDataset(os.path.join(self.filename, "*.nc"), 'r')
        # elif self.filename.endswith("*"):
        #     self.ncfile = MFDataset(self.filename + ".nc", 'r')
        # else:
        #     raise ValueError("Unable to decode file " + str(self.filename))

    def find_time_and_step(self,index_t):
        try:
            if "time" in self.ds_instant.variables and "time" in self.ds_accum.variables and "step" in self.ds_accum.variables:
                array = np.asarray(self.ds_accum.variables["time"].data)
                idx = np.where(array <= self.ds_instant.variables["time"].data[int(index_t)])

                if len(idx[0]) == 0:
                    raise ValueError("Time "+self.ds_instant.variables["time"].data[int(index_t)]+" was not found")

                nearest_t_index = (np.abs(array[idx] - self.ds_instant.variables["time"].data[int(index_t)])).argmin()

                td = datetime.utcfromtimestamp(
                    self.ds_instant.variables["time"].data[int(index_t)]) - datetime.utcfromtimestamp(
                    self.ds_accum.variables["time"].data[int(nearest_t_index)])

                array = np.asarray(self.ds_accum.variables["step"].data)
                nearest_step_index = (np.abs(array - divmod(td.seconds, 3600)[0])).argmin()

                #logging.debug("On cherche :",datetime.utcfromtimestamp(self.ds_instep.variables["time"].data[int(index_t)]) )
                #logging.debug("ref Time",datetime.utcfromtimestamp(self.ds_step.variables["time"].data[int(nearest_t_index)]))
                step = timedelta(hours=self.ds_accum.variables["step"].data[int(nearest_step_index)])
                #logging.debug("On trouve",datetime.utcfromtimestamp(self.ds_step.variables["time"].data[int(nearest_t_index)])+step)

                return nearest_t_index, nearest_step_index
            else:
                logging.debug("No variables found for '" + str(VariableDefinition.LONG_NAME['time']) + "'")
                raise (VariableNameError("ECMWFReader",
                                         "No variables found for '" + str(
                                             VariableDefinition.LONG_NAME['time']) + "'",
                                         1000))
        except Exception as ex:
            logging.debug("Error '" + str(ex) + "'")
            raise (VariableNameError("ECMWFReader", "An error occured : '" + str(ex) + "'", 1000))

    def close(self):
        self.ds_instant.close()

    def is_regular_grid(self):
        return True

    def get_x_size(self):
        try:
            if "longitude" in self.ds_instant.variables:
                return self.ds_instant.dimensions["longitude"]
            else:
                logging.debug("No variables found for '" + str(VariableDefinition.LONG_NAME['longitude']) + "'")
                raise (VariableNameError("ECMWFReader",
                                         "No variables found for '" + str(
                                             VariableDefinition.LONG_NAME['longitude']) + "'",
                                         1000))
        except Exception as ex:
            logging.debug("Error '" + str(ex) + "'")
            raise (VariableNameError("ECMWFReader", "An error occured : '" + str(ex) + "'", 1000))

    def get_y_size(self):
        try:
            if "latitude" in self.ds_instant.variables:
                return self.ds_instant.dimensions["latitude"]
            else:
                logging.debug("No variables found for '" + str(VariableDefinition.LONG_NAME['latitude']) + "'")
                raise (VariableNameError("ECMWFReader",
                                         "No variables found for '" + str(
                                             VariableDefinition.LONG_NAME['latitude']) + "'",
                                         1000))
        except Exception as ex:
            logging.debug("Error '" + str(ex) + "'")
            raise (VariableNameError("ECMWFReader", "An error occured : '" + str(ex) + "'", 1000))

    def get_t_size(self):
        try:
            if "time" in self.ds_instant.variables:
                return self.ds_instant.dimensions["time"]
            else:
                logging.debug("No variables found for '" + str(VariableDefinition.LONG_NAME['time']) + "'")
                raise (VariableNameError("ECMWFReader",
                                         "No variables found for '" + str(
                                             VariableDefinition.LONG_NAME['time']) + "'",
                                         1000))
        except Exception as ex:
            logging.debug("Error '" + str(ex) + "'")
            raise (VariableNameError("ECMWFReader", "An error occured : '" + str(ex) + "'", 1000))

    def read_axis_x(self,xmin,xmax,ymin,ymax):
        return np.array(sorted(self.new_lon)[xmin:xmax])

    def read_axis_y(self,xmin,xmax,ymin,ymax):
        try:
            if "latitude" in self.ds_instant.variables:
                return self.ds_instant.variables['latitude'].data[ymin:ymax]
            else:
                logging.debug("No variables found for '" + str(VariableDefinition.LONG_NAME['latitude']) + "'")
                raise (VariableNameError("ECMWFReader",
                                         "No variables found for '" + str(
                                             VariableDefinition.LONG_NAME['latitude']) + "'",
                                         1000))
        except Exception as ex:
            logging.debug("Error '" + str(ex) + "'")
            raise (VariableNameError("ECMWFReader", "An error occured : '" + str(ex) + "'", 1000))

    def read_axis_t(self,tmin,tmax,timestamp):
        try:
            if "time" in self.ds_instant.variables:
                data = self.ds_instant.variables['time'].data[tmin:tmax]
                result = [datetime.utcfromtimestamp(t) \
                          for t in data]

                if timestamp == 1:
                    return [(t - TimeCoverage.TIME_DATUM).total_seconds() \
                            for t in result];
                else:
                    return result;
            else:
                logging.debug("No variables found for '" + str(VariableDefinition.LONG_NAME['time']) + "'")
                raise (VariableNameError("ECMWFReader",
                                         "No variables found for '" + str(
                                             VariableDefinition.LONG_NAME['time']) + "'",
                                         1000))
        except Exception as ex:
            logging.debug("Error '" + str(ex) + "'")
            raise (VariableNameError("ECMWFReader", "An error occured : '" + str(ex) + "'", 1000))

    # Variables
    def read_variable_longitude(self,xmin,xmax,ymin,ymax):
        return self.read_axis_x()

    def read_variable_latitude(self,xmin,xmax,ymin,ymax):
        return self.read_axis_y()

    def read_variable_time(self,tmin,tmax):
        return self.read_axis_t(tmin,tmax,timestamp=0)

    def read_variable_2D_land_binary_mask(self,xmin,xmax,ymin,ymax):
        try:
            if "lsm" in self.ds_instant.variables:
                mask = np.ma.filled(self.ds_instant.variables["lsm"].data[0, ymin:ymax, :], fill_value=np.nan)
                return np.take_along_axis(mask, self.new_indexes[ymin:ymax], axis=1)[:, xmin:xmax]
            else:
                logging.debug("No variables found for '" + str(VariableDefinition.LONG_NAME['2d_land_binary_mask']) + "'")
                raise (VariableNameError("ECMWFReader",
                                         "No variables found for '" + str(
                                             VariableDefinition.LONG_NAME['2d_land_binary_mask']) + "'",
                                         1000))
        except Exception as ex:
            logging.debug("Error '" + str(ex) + "'")
            raise (VariableNameError("ECMWFReader", "An error occured : '" + str(ex) + "'", 1000))

    def read_variable_2D_land_binary_mask_at_time(self,index_t, xmin, xmax, ymin, ymax):
        try:
            if "lsm" in self.ds_instant.variables:
                mask = np.ma.filled(self.ds_instant.variables["lsm"].data[int(index_t), ymin:ymax, :],
                                    fill_value=np.nan)
                return np.take_along_axis(mask, self.new_indexes[ymin:ymax], axis=1)[:, xmin:xmax]
            else:
                logging.debug(
                    "No variables found for '" + str(VariableDefinition.LONG_NAME['2d_land_binary_mask']) + "'")
                raise (VariableNameError("ECMWFReader",
                                         "No variables found for '" + str(
                                             VariableDefinition.LONG_NAME['2d_land_binary_mask']) + "'",
                                         1000))
        except Exception as ex:
            logging.debug("Error '" + str(ex) + "'")
            raise (VariableNameError("ECMWFReader", "An error occured : '" + str(ex) + "'", 1000))

    def read_variable_2D_sea_binary_mask(self,xmin,xmax,ymin,ymax):
        try:
            if "lsm" in self.ds_instant.variables:
                mask = np.ma.filled(self.ds_instant.variables["lsm"].data[0, ymin:ymax, :], fill_value=np.nan)
                mask += 1.0  # inverse le mask
                mask %= 2  # inverse le mask
                return np.take_along_axis(mask, self.new_indexes[ymin:ymax], axis=1)[:, xmin:xmax]
            else:
                logging.debug(
                    "No variables found for '" + str(VariableDefinition.LONG_NAME['2d_sea_binary_mask']) + "'")
                raise (VariableNameError("ECMWFReader",
                                         "No variables found for '" + str(
                                             VariableDefinition.LONG_NAME['2d_sea_binary_mask']) + "'",
                                         1000))
        except Exception as ex:
            logging.debug("Error '" + str(ex) + "'")
            raise (VariableNameError("ECMWFReader", "An error occured : '" + str(ex) + "'", 1000))

    def read_variable_3D_sea_binary_mask_at_time(self, index_t,xmin,xmax,ymin,ymax):
        try:
            if "lsm" in self.ds_instant.variables:
                mask = np.ma.filled(self.ds_instant.variables["lsm"].data[int(index_t), int(ymin):int(ymax), :],
                                    fill_value=np.nan)
                mask += 1.0  # inverse le mask
                mask %= 2  # inverse le mask
                return np.take_along_axis(mask, self.new_indexes[ymin:ymax], axis=1)[:, xmin:xmax]
            else:
                logging.debug(
                    "No variables found for '" + str(VariableDefinition.LONG_NAME['sea_binary_mask']) + "'")
                raise (VariableNameError("ECMWFReader",
                                         "No variables found for '" + str(
                                             VariableDefinition.LONG_NAME['sea_binary_mask']) + "'",
                                         1000))
        except Exception as ex:
            logging.debug("Error '" + str(ex) + "'")
            raise (VariableNameError("ECMWFReader", "An error occured : '" + str(ex) + "'", 1000))

    def read_variable_3D_land_binary_mask_at_time(self, index_t,xmin,xmax,ymin,ymax):
        try:
            if "lsm" in self.ds_instant.variables:
                mask = np.ma.filled(self.ds_instant.variables["lsm"].data[int(index_t), int(ymin):int(ymax), :],
                                    fill_value=np.nan)
                return np.take_along_axis(mask, self.new_indexes[ymin:ymax], axis=1)[:, xmin:xmax]
            else:
                logging.debug(
                    "No variables found for '" + str(VariableDefinition.LONG_NAME['land_binary_mask']) + "'")
                raise (VariableNameError("ECMWFReader",
                                         "No variables found for '" + str(
                                             VariableDefinition.LONG_NAME['land_binary_mask']) + "'",
                                         1000))
        except Exception as ex:
            logging.debug("Error '" + str(ex) + "'")
            raise (VariableNameError("ECMWFReader", "An error occured : '" + str(ex) + "'", 1000))

    #################
    # METEO
    # 2D
    #################

    def read_variable_rainfall_amount_at_time(self, index_t,xmin,xmax,ymin,ymax):
        try:
            if "tp" in self.ds_accum.variables:
                nearest_t_index, nearest_step_index = self.find_time_and_step(index_t)
                data = np.ma.filled(self.ds_accum.variables["tp"].data[int(nearest_t_index), int(nearest_step_index),
                                    int(ymin):int(ymax), :], fill_value=np.nan)
                return np.take_along_axis(data, self.new_indexes[ymin:ymax], axis=1)[:, xmin:xmax]
            else:
                logging.debug(
                    "No variables found for '" + str(VariableDefinition.LONG_NAME['rainfall_amount']) + "'")
                raise (VariableNameError("ECMWFReader",
                                         "No variables found for '" + str(
                                             VariableDefinition.LONG_NAME['rainfall_amount']) + "'",
                                         1000))
        except Exception as ex:
            logging.debug("Error '" + str(ex) + "'")
            raise (VariableNameError("ECMWFReader", "An error occured : '" + str(ex) + "'", 1000))

    #################
    # METEO
    # Sea surface
    #################

    def read_variable_surface_air_pressure_at_time(self, index_t,xmin,xmax,ymin,ymax):
        try:
            if "sp" in self.ds_instant.variables:
                data = np.ma.filled(self.ds_instant.variables["sp"].data[int(index_t), int(ymin):int(ymax), :],
                                    fill_value=np.nan)
                data *= 0.01  # Pa to hPa
                return np.take_along_axis(data, self.new_indexes[ymin:ymax], axis=1)[:, xmin:xmax]
            else:
                logging.debug(
                    "No variables found for '" + str(VariableDefinition.LONG_NAME['surface_air_pressure']) + "'")
                raise (VariableNameError("ECMWFReader",
                                         "No variables found for '" + str(
                                             VariableDefinition.LONG_NAME['surface_air_pressure']) + "'",
                                         1000))
        except Exception as ex:
            logging.debug("Error '" + str(ex) + "'")
            raise (VariableNameError("ECMWFReader", "An error occured : '" + str(ex) + "'", 1000))

    def read_variable_sea_surface_air_pressure_at_time(self,index_t,xmin,xmax,ymin,ymax):
        try:
            if "msl" in self.ds_instant.variables:
                data = np.ma.filled(self.ds_instant.variables["msl"].data[int(index_t), int(ymin):int(ymax), :],
                                    fill_value=np.nan)
                data *= 0.01  # Pa to hPa
                return np.take_along_axis(data, self.new_indexes[ymin:ymax], axis=1)[:, xmin:xmax]
            else:
                logging.debug(
                    "No variables found for '" + str(VariableDefinition.LONG_NAME['sea_surface_air_pressure']) + "'")
                raise (VariableNameError("ECMWFReader",
                                         "No variables found for '" + str(
                                             VariableDefinition.LONG_NAME['sea_surface_air_pressure']) + "'",
                                         1000))
        except Exception as ex:
            logging.debug("Error '" + str(ex) + "'")
            raise (VariableNameError("ECMWFReader", "An error occured : '" + str(ex) + "'", 1000))

    def read_variable_surface_downward_sensible_heat_flux_at_time(self, index_t,xmin,xmax,ymin,ymax):
        try:
            if "sshf" in self.ds_accum.variables:
                nearest_t_index, nearest_step_index = self.find_time_and_step(index_t)
                data = np.ma.filled(self.ds_accum.variables["sshf"].data[int(nearest_t_index), int(nearest_step_index),
                                    int(ymin):int(ymax), :], fill_value=np.nan)
                return np.take_along_axis(data, self.new_indexes[ymin:ymax], axis=1)[:, xmin:xmax] #/ 86400.00
            else:
                logging.debug(
                    "No variables found for '" + str(VariableDefinition.LONG_NAME['surface_downward_sensible_heat_flux']) + "'")
                raise (VariableNameError("ECMWFReader",
                                         "No variables found for '" + str(
                                             VariableDefinition.LONG_NAME['surface_downward_sensible_heat_flux']) + "'",
                                         1000))
        except Exception as ex:
            logging.debug("Error '" + str(ex) + "'")
            raise (VariableNameError("ECMWFReader", "An error occured : '" + str(ex) + "'", 1000))

    def read_variable_surface_downward_latent_heat_flux_at_time(self, index_t,xmin,xmax,ymin,ymax):
        try:
            if "slhf" in self.ds_accum.variables:
                nearest_t_index, nearest_step_index = self.find_time_and_step(index_t)
                data = np.ma.filled(self.ds_accum.variables["slhf"].data[int(nearest_t_index), int(nearest_step_index),
                                    int(ymin):int(ymax), :], fill_value=np.nan)
                return np.take_along_axis(data, self.new_indexes[ymin:ymax], axis=1)[:, xmin:xmax] #/ 86400.00
            else:
                logging.debug(
                    "No variables found for '" + str(
                        VariableDefinition.LONG_NAME['surface_downward_latent_heat_flux']) + "'")
                raise (VariableNameError("ECMWFReader",
                                         "No variables found for '" + str(
                                             VariableDefinition.LONG_NAME['surface_downward_latent_heat_flux']) + "'",
                                         1000))
        except Exception as ex:
            logging.debug("Error '" + str(ex) + "'")
            raise (VariableNameError("ECMWFReader", "An error occured : '" + str(ex) + "'", 1000))

    def read_variable_surface_air_temperature_at_time(self, index_t,xmin,xmax,ymin,ymax):
        try:
            if "t2m" in self.ds_instant.variables:
                data = np.ma.filled(self.ds_instant.variables["t2m"].data[int(index_t), int(ymin):int(ymax), :] - 273.15,
                                    fill_value=np.nan)
                return np.take_along_axis(data, self.new_indexes[ymin:ymax], axis=1)[:, xmin:xmax]
            else:
                logging.debug(
                    "No variables found for '" + str(
                        VariableDefinition.LONG_NAME['surface_air_temperature']) + "'")
                raise (VariableNameError("ECMWFReader",
                                         "No variables found for '" + str(
                                             VariableDefinition.LONG_NAME['surface_air_temperature']) + "'",
                                         1000))
        except Exception as ex:
            logging.debug("Error '" + str(ex) + "'")
            raise (VariableNameError("ECMWFReader", "An error occured : '" + str(ex) + "'", 1000))

    def read_variable_dew_point_temperature_at_time(self, index_t,xmin,xmax,ymin,ymax):
        try:
            if "d2m" in self.ds_instant.variables:
                data = np.ma.filled(self.ds_instant.variables["d2m"].data[int(index_t), int(ymin):int(ymax), :] - 273.15,
                                    fill_value=np.nan)
                return np.take_along_axis(data, self.new_indexes[ymin:ymax], axis=1)[:, xmin:xmax]
            else:
                logging.debug(
                    "No variables found for '" + str(
                        VariableDefinition.LONG_NAME['dew_point_temperature']) + "'")
                raise (VariableNameError("ECMWFReader",
                                         "No variables found for '" + str(
                                             VariableDefinition.LONG_NAME['dew_point_temperature']) + "'",
                                         1000))
        except Exception as ex:
            logging.debug("Error '" + str(ex) + "'")
            raise (VariableNameError("ECMWFReader", "An error occured : '" + str(ex) + "'", 1000))

    def read_variable_surface_downward_solar_radiation_at_time(self, index_t, xmin, xmax, ymin, ymax):
        try:
            if "ssrd" in self.ds_accum.variables:
                nearest_t_index, nearest_step_index = self.find_time_and_step(index_t)
                data = np.ma.filled(self.ds_accum.variables["ssrd"].data[int(nearest_t_index), int(nearest_step_index),
                                    int(ymin):int(ymax), :], fill_value=np.nan)
                return np.take_along_axis(data, self.new_indexes[ymin:ymax], axis=1)[:, xmin:xmax] #/ 86400.00
            else:
                logging.debug(
                    "No variables found for '" + str(
                        VariableDefinition.LONG_NAME['surface_downwards_solar_radiation']) + "'")
                raise (VariableNameError("ECMWFReader",
                                         "No variables found for '" + str(
                                             VariableDefinition.LONG_NAME['surface_downwards_solar_radiation']) + "'",
                                         1000))
        except Exception as ex:
            logging.debug("Error '" + str(ex) + "'")
            raise (VariableNameError("ECMWFReader", "An error occured : '" + str(ex) + "'", 1000))

    def read_variable_surface_downward_thermal_radiation_at_time(self, index_t, xmin, xmax, ymin, ymax):
        try:
            if "strd" in self.ds_accum.variables:
                nearest_t_index, nearest_step_index = self.find_time_and_step(index_t)
                data = np.ma.filled(self.ds_accum.variables["strd"].data[int(nearest_t_index), int(nearest_step_index),
                                    int(ymin):int(ymax), :], fill_value=np.nan)
                return np.take_along_axis(data, self.new_indexes[ymin:ymax], axis=1)[:, xmin:xmax] #/ 86400.00
            else:
                logging.debug(
                    "No variables found for '" + str(
                        VariableDefinition.LONG_NAME['surface_downwards_thermal_radiation']) + "'")
                raise (VariableNameError("ECMWFReader",
                                         "No variables found for '" + str(
                                             VariableDefinition.LONG_NAME['surface_downwards_thermal_radiation']) + "'",
                                         1000))
        except Exception as ex:
            logging.debug("Error '" + str(ex) + "'")
            raise (VariableNameError("ECMWFReader", "An error occured : '" + str(ex) + "'", 1000))

    def read_variable_surface_solar_radiation_at_time(self, index_t,xmin,xmax,ymin,ymax):
        try:
            if "ssr" in self.ds_accum.variables:
                nearest_t_index, nearest_step_index = self.find_time_and_step(index_t)
                data = np.ma.filled(self.ds_accum.variables["ssr"].data[int(nearest_t_index), int(nearest_step_index),
                                    int(ymin):int(ymax), :], fill_value=np.nan)
                return np.take_along_axis(data, self.new_indexes[ymin:ymax], axis=1)[:, xmin:xmax] #/ 86400.00
            else:
                logging.debug(
                    "No variables found for '" + str(
                        VariableDefinition.LONG_NAME['surface_solar_radiation']) + "'")
                raise (VariableNameError("ECMWFReader",
                                         "No variables found for '" + str(
                                             VariableDefinition.LONG_NAME['surface_solar_radiation']) + "'",
                                         1000))
        except Exception as ex:
            logging.debug("Error '" + str(ex) + "'")
            raise (VariableNameError("ECMWFReader", "An error occured : '" + str(ex) + "'", 1000))

    def read_variable_surface_thermal_radiation_at_time(self, index_t,xmin,xmax,ymin,ymax):
        try:
            if "str" in self.ds_accum.variables:
                nearest_t_index, nearest_step_index = self.find_time_and_step(index_t)
                data = np.ma.filled(self.ds_accum.variables["str"].data[int(nearest_t_index), int(nearest_step_index),
                                    int(ymin):int(ymax), :], fill_value=np.nan)
                return np.take_along_axis(data, self.new_indexes[ymin:ymax], axis=1)[:, xmin:xmax] #/ 86400.00
            else:
                logging.debug(
                    "No variables found for '" + str(
                        VariableDefinition.LONG_NAME['surface_thermal_radiation']) + "'")
                raise (VariableNameError("ECMWFReader",
                                         "No variables found for '" + str(
                                             VariableDefinition.LONG_NAME['surface_thermal_radiation']) + "'",
                                         1000))
        except Exception as ex:
            logging.debug("Error '" + str(ex) + "'")
            raise (VariableNameError("ECMWFReader", "An error occured : '" + str(ex) + "'", 1000))

    #################
    # METEO
    # At 10 m
    #################

    def read_variable_wind_10m_at_time(self, index_t,xmin,xmax,ymin,ymax):
        try:
            if "u10" in self.ds_instant.variables and "v10" in self.ds_instant.variables:
                u = np.ma.filled(self.ds_instant.variables["u10"].data[int(index_t), int(ymin):int(ymax), :],
                                 fill_value=np.nan)
                v = np.ma.filled(self.ds_instant.variables["v10"].data[int(index_t), int(ymin):int(ymax), :],
                                 fill_value=np.nan)
                new_u = np.take_along_axis(u, self.new_indexes[ymin:ymax], axis=1)
                new_v = np.take_along_axis(v, self.new_indexes[ymin:ymax], axis=1)
                return [new_u[:, xmin:xmax], new_v[:, xmin:xmax]]
            else:
                logging.debug("No variables found for \'Wind 10m\'")
                raise (VariableNameError("ECMWFReader",
                                         "No variables found for \'Wind 10m\'",
                                         1000))
        except Exception as ex:
            logging.debug("Error '" + str(ex) + "'")
            raise (VariableNameError("ECMWFReader", "An error occured : '" + str(ex) + "'", 1000))


