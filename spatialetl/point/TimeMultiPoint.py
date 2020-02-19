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
from spatialetl.point.MultiPoint import MultiPoint
from spatialetl.operator.interpolator.InterpolatorCore import time_1d_interpolation
from numpy import int,int32,int64
import numpy as np
from datetime import datetime,timedelta,timezone
from spatialetl.utils.logger import logging
import pandas
import math
import cftime

class TimeMultiPoint(MultiPoint):
    """"""

    TIME_DATUM = datetime(1970, 1, 1)
    TIME_DELTA = timedelta(minutes=5)
    TIME_INTERPOLATION_METHOD = "linear"

    def __init__(self,myReader,start_time=None,end_time=None,freq=None,time_range=None):
        MultiPoint.__init__(self, myReader)

        self.source_axis_t = self.reader.read_axis_t(timestamp=0)
        self.source_t_size = np.shape(self.source_axis_t)[0]
        tmin = 0
        tmax = self.source_t_size
        zero_delta = timedelta(minutes=00)

        self.temporal_resampling = False

        if time_range is not None:
            self.target_axis_t = time_range
            self.target_t_size = np.shape(self.target_axis_t)[0]
            self.temporal_resampling = True

        if start_time is not None:

            if type(start_time) == datetime or type(start_time) == cftime._cftime.real_datetime:
                time = start_time
            elif type(start_time) == str:
                try:
                    time = datetime.strptime(start_time, '%Y-%m-%d %H:%M:%S')
                except ValueError as ex:
                    raise ValueError("start_time is not well formated " + str(ex))
            else:
                raise ValueError("start_time have to be string or datetime. Found " + str(type(start_time)))

            nearest_t_index = (np.abs(np.asarray(self.source_axis_t) - time)).argmin()

            if time - self.source_axis_t[nearest_t_index] == zero_delta or abs(
                    time - self.source_axis_t[nearest_t_index]) < TimeMultiPoint.TIME_DELTA:
                tmin = nearest_t_index
            else:
                raise ValueError(str(time) + " not found. Maybe the TimeMultiPoint.TIME_DELTA (" + str(
                    TimeMultiPoint.TIME_DELTA) + ") is too small or the date is out the range.")

        if end_time is not None:

            if type(end_time) == datetime or type(end_time) == cftime._cftime.real_datetime:
                time = end_time
            elif type(end_time) == str:
                try:
                    time = datetime.strptime(end_time, '%Y-%m-%d %H:%M:%S')
                except ValueError as ex:
                    raise ValueError("end_time is not well formated " + str(ex))
            else:
                raise ValueError("end_time have to be string or datetime. Found " + str(type(end_time)))

            nearest_t_index = (np.abs(np.asarray(self.source_axis_t) - time)).argmin()

            if time - self.source_axis_t[nearest_t_index] == zero_delta or abs(
                    time - self.source_axis_t[nearest_t_index]) < TimeMultiPoint.TIME_DELTA:
                tmax = nearest_t_index + 1
            else:
                raise ValueError(str(time) + " not found. Maybe the TimeMultiPoint.TIME_DELTA (" + str(
                    TimeMultiPoint.TIME_DELTA) + ") is too small or the date is out the range.")

        if freq is not None:
            self.temporal_resampling = True
            self.target_axis_t = pandas.date_range(start=self.source_axis_t[tmin],
                                                          end=self.source_axis_t[tmax - 1],
                                                          freq=freq).to_pydatetime();
            self.target_t_size = np.shape(self.target_axis_t)[0]
        else:
            self.target_axis_t = self.source_axis_t[tmin:tmax]
            self.target_t_size = tmax - tmin

    # Axis
    def read_axis_t(self,timestamp=0,type="target"):

        if  type=="source" and timestamp==1:
            return self.reader.read_axis_t(timestamp)
        elif type=="source":
            return self.source_axis_t
        else:
            return self.target_axis_t

    def get_t_size(self,type="target"):

        if type=="source":
            return self.source_t_size
        else:
            return self.target_t_size

    def find_time_index(self, t):
        """Retourne l'index de la date la plus proche à TIME_DELTA_MIN prêt.
    @type t: datetime ou int
    @param t: date souhaitée ou l'index de la date souhaitée
    @return:  l'index de la date la plus proche à TIME_DELTA_MIN prêt ou une erreur si aucune date n'a pu être trouvée."""

        indexes_t = []
        zero_delta = timedelta(seconds=0)

        if type(t) == int or type(t) == int32 or type(t)== int64:

            if t < 0 or t >= self.get_t_size(type="source"):
                raise ValueError("Time index have to range between 0 and " + str(
                    self.get_t_size(type="source") - 1) + ". Actually Time index = " + str(t))

            indexes_t.append(int(t));

        elif type(t) == datetime or type(t) == cftime._cftime.real_datetime:

            logging.debug("[TimeMultiPoint][find_time_index()] Looking for : "+str(t))

            array = np.asarray(self.read_axis_t(timestamp=1,type="source"))
            #X = np.abs(self.read_axis_t(timestamp=1,raw=1) - t.replace(tzinfo=timezone.utc).timestamp())

            X = np.abs(array - t)
            #index_t = (np.abs(array - t.replace(tzinfo=timezone.utc).timestamp())).argmin()
            idx = np.where(X <= TimeMultiPoint.TIME_DELTA)

            for index in range(np.shape(idx)[1]):
                index_t = idx[0][index]
                indexes_t.append(int(index_t))
                logging.debug("[TimeMultiPoint][find_time_index()] Found : " + str(self.source_axis_t[index_t]))

            if not indexes_t:
                raise ValueError("" + str(t) + " was not found. Maybe the TimeMultiPoint.TIME_DELTA (" + str(
                    TimeMultiPoint.TIME_DELTA) + ") is too small or the date is out the range.")

        else:
            raise ValueError("" + str(t) + " have to be an integer or a datetime.")

        logging.debug("[TimeMultiPoint][find_time_index()] Found " + str(len(indexes_t)) + " candidate datetime(s)")

        # On retourne le tableau d'index
        return np.array(np.unique(indexes_t))

    def interpolate_all_times(self, values):

        source_index = pandas.DatetimeIndex(self.read_axis_t(type="source"))
        target_index = pandas.DatetimeIndex(self.read_axis_t(type="target"))

        for index_x in range(0,self.get_nb_points()):

            data = pandas.DataFrame({'point_'+str(index_x): pandas.Series(values[index_x], index=source_index)})

        # we process time record (drop duplicate...)
        duplicates = np.where(data.index.duplicated() == True)[0]

        count = np.shape(duplicates)[0]
        if count > 0:
            logging.warning(
                '[TimeMultiPoint] ' + str(count) + ' dates are duplicated. We drop them by keeping the first.')
            logging.debug("[TimeMultiPoint] current file= " + str(self.reader.filename))
            logging.debug("[TimeMultiPoint] duplicates index=" + str(duplicates))
            data = data.groupby(data.index).first()

        # we sort by dates
        data = data.sort_index()

        if TimeMultiPoint.TIME_INTERPOLATION_METHOD == None:
            data = data.reindex(target_index, fill_value=np.nan)

        if TimeMultiPoint.TIME_INTERPOLATION_METHOD == "linear":
            # Interpolate
            tmp = data.reindex(target_index, fill_value=np.nan)
            data = tmp.interpolate(method=TimeMultiPoint.TIME_INTERPOLATION_METHOD, axis=0)
        else:
            # Fill gap
            data = data.reindex(target_index, method=TimeMultiPoint.TIME_INTERPOLATION_METHOD, fill_value=np.nan,
                                tolerance=TimeMultiPoint.TIME_DELTA);

        return data.as_matrix()

    def interpolate_time(self,date,indexes_t,layers):

        results = np.zeros([self.get_nb_points()])
        results[:] = np.NAN

        targetTime = [date.replace(tzinfo=timezone.utc).timestamp()]
        rawTime = self.read_axis_t(timestamp=1, type="source")

        candidateTimes = np.zeros([len(indexes_t)])
        candidateValues = np.zeros([len(indexes_t)])

        for x in range(0,self.get_nb_points()):

            candidateValues[:] = np.nan
            candidateTimes[:] = np.nan

            for t in range(0, len(indexes_t)):
                candidateTimes[t] = datetime.timestamp(rawTime[indexes_t[t]])
                candidateValues[t] = layers[t][x]

            # We remove NaN values
            nan_indexes = np.isnan(candidateValues)
            finalCandidateTimes = candidateTimes[~nan_indexes]
            finalCandidateValues = candidateValues[~nan_indexes]

            if len(finalCandidateValues) > 1:
                # We have more than 1 value, we interpole them
                results[x] = time_1d_interpolation(finalCandidateTimes, targetTime, finalCandidateValues,TimeMultiPoint.TIME_INTERPOLATION_METHOD)
            else:
                # We take the only value
                results[x] = candidateValues[0]

        return results

    # Scalar
    def read_variable_longitude_at_time(self, date):
        index_t = self.find_time_index(date)

        if len(index_t) > 1:
            layers = np.zeros([len(index_t), self.get_nb_points()])
            layers[::] = np.NAN

            for t in range(0, len(index_t)):
                layers[t] = self.reader.read_variable_longitude_at_time(index_t[t])

            data = self.interpolate_time(date,index_t,layers)

        else:
            data = self.reader.read_variable_longitude_at_time(index_t[0])

        return data

    def read_variable_latitude_at_time(self, date):
        index_t = self.find_time_index(date)

        if len(index_t) > 1:
            layers = np.zeros([len(index_t), self.get_nb_points()])
            layers[::] = np.NAN

            for t in range(0, len(index_t)):
                layers[t] = self.reader.read_variable_latitude_at_time(index_t[t])

            data = self.interpolate_time(date,index_t,layers)

        else:
            data = self.reader.read_variable_latitude_at_time(index_t[0])

        return data

    # HYDRO
    def read_variable_sea_surface_temperature_at_time(self, date):
        index_t = self.find_time_index(date)

        if len(index_t) > 1:
            layers = np.zeros([len(index_t), self.get_nb_points()])
            layers[::] = np.NAN

            for t in range(0, len(index_t)):
                layers[t] = self.reader.read_variable_sea_surface_temperature_at_time(index_t[t])

            data = self.interpolate_time(date,index_t,layers)

        else:
            data = self.reader.read_variable_sea_surface_temperature_at_time(index_t[0])

        return data

    def read_variable_sea_water_temperature_at_ground_level_at_time(self, date):
        index_t = self.find_time_index(date)

        if len(index_t) > 1:
            layers = np.zeros([len(index_t), self.get_nb_points()])
            layers[::] = np.NAN

            for t in range(0, len(index_t)):
                layers[t] = self.reader.read_variable_sea_water_temperature_at_ground_level_at_time(index_t[t])

            data = self.interpolate_time(date, index_t, layers)

        else:
            data = self.reader.read_variable_sea_water_temperature_at_ground_level_at_time(index_t[0])

        return data

    def read_variable_sea_surface_salinity_at_time(self, date):
        index_t = self.find_time_index(date)

        if len(index_t) > 1:
            layers = np.zeros([len(index_t), self.get_nb_points()])
            layers[::] = np.NAN

            for t in range(0, len(index_t)):
                layers[t] = self.reader.read_variable_sea_surface_salinity_at_time(index_t[t])

            data = self.interpolate_time(date,index_t,layers)

        else:
            data = self.reader.read_variable_sea_surface_salinity_at_time(index_t[0])

        return data

    def read_variable_sea_water_salinity_at_ground_level_at_time(self, date):
        index_t = self.find_time_index(date)

        if len(index_t) > 1:
            layers = np.zeros([len(index_t), self.get_nb_points()])
            layers[::] = np.NAN

            for t in range(0, len(index_t)):
                layers[t] = self.reader.read_variable_sea_water_salinity_at_ground_level_at_time(index_t[t])

            data = self.interpolate_time(date, index_t, layers)

        else:
            data = self.reader.read_variable_sea_water_salinity_at_ground_level_at_time(index_t[0])

        return data

    def read_variable_sea_water_pressure_at_sea_water_surface_at_time(self, date):
        index_t = self.find_time_index(date)

        if len(index_t) > 1:
            layers = np.zeros([len(index_t), self.get_nb_points()])
            layers[::] = np.NAN

            for t in range(0, len(index_t)):
                layers[t] = self.reader.read_variable_sea_water_pressure_at_sea_water_surface_at_time(index_t[t])

            data = self.interpolate_time(date,index_t,layers)

        else:
            data = self.reader.read_variable_sea_water_pressure_at_sea_water_surface_at_time(index_t[0])

        return data

    def read_variable_sea_surface_height_above_mean_sea_level_at_time(self, date):
        index_t = self.find_time_index(date)

        if len(index_t) > 1:
            layers = np.zeros([len(index_t), self.get_nb_points()])
            layers[::] = np.NAN

            for t in range(0, len(index_t)):
                layers[t] = self.reader.read_variable_sea_surface_height_above_mean_sea_level_at_time(index_t[t])

            data = self.interpolate_time(date, index_t, layers)

        else:
            data = self.reader.read_variable_sea_surface_height_above_mean_sea_level_at_time(index_t[0])

        return data

    def read_variable_sea_surface_height_above_geoid_at_time(self, date):
        index_t = self.find_time_index(date)

        if len(index_t) > 1:
            layers = np.zeros([len(index_t), self.get_nb_points()])
            layers[::] = np.NAN

            for t in range(0, len(index_t)):
                layers[t] = self.reader.read_variable_sea_surface_height_above_geoid_at_time(index_t[t])

            data = self.interpolate_time(date, index_t, layers)

        else:
            data = self.reader.read_variable_sea_surface_height_above_geoid_at_time(index_t[0])

        return data

    def read_variable_sea_surface_density_at_time(self, date):
        index_t = self.find_time_index(date)

        if len(index_t) > 1:
            layers = np.zeros([len(index_t), self.get_nb_points()])
            layers[::] = np.NAN

            for t in range(0, len(index_t)):
                layers[t] = self.reader.read_variable_sea_surface_density_at_time(index_t[t])

            data = self.interpolate_time(date,index_t,layers)

        else:
            data = self.reader.read_variable_sea_surface_density_at_time(index_t[0])

        return data

    def read_variable_sea_water_turbidity_at_time(self, date):
        index_t = self.find_time_index(date)

        if len(index_t) > 1:
            layers = np.zeros([len(index_t), self.get_nb_points()])
            layers[::] = np.NAN

            for t in range(0, len(index_t)):
                layers[t] = self.reader.read_variable_sea_water_turbidity_at_time(index_t[t])

            data = self.interpolate_time(date,index_t,layers)

        else:
            data = self.reader.read_variable_sea_water_turbidity_at_time(index_t[0])

        return data

    def read_variable_sea_water_electrical_conductivity_at_time(self, date):
        index_t = self.find_time_index(date)

        if len(index_t) > 1:
            layers = np.zeros([len(index_t), self.get_nb_points()])
            layers[::] = np.NAN

            for t in range(0, len(index_t)):
                layers[t] = self.reader.read_variable_sea_water_electrical_conductivity_at_time(index_t[t])

            data = self.interpolate_time(date,index_t,layers)

        else:
            data = self.reader.read_variable_sea_water_electrical_conductivity_at_time(index_t[0])

        return data


    def read_variable_barotropic_sea_water_velocity_at_time(self, date):
        index_t = self.find_time_index(date)

        data = np.zeros([2, self.get_nb_points()])
        data[::] = np.NAN

        if len(index_t) > 1:
            layers = np.zeros([len(index_t), 2, self.get_nb_points()])
            layers[::] = np.NAN

            for t in range(0, len(index_t)):
                comp = self.reader.read_variable_barotropic_sea_water_velocity_at_time(index_t[t])
                layers[t][0] = comp[0]
                layers[t][1] = comp[1]

            data[0] = self.interpolate_time(date, index_t, layers[:, 0, :])
            data[1] = self.interpolate_time(date, index_t, layers[:, 1, :])

        else:
            data = self.reader.read_variable_barotropic_sea_water_velocity_at_time(index_t[0])

        return data

    def read_variable_barotropic_sea_water_speed_at_time(self, date):
        comp = self.read_variable_barotropic_sea_water_velocity_at_time(date)

        data = np.zeros([self.get_nb_points()])
        data[::] = np.NAN

        for index_x in range(0, self.get_nb_points()):
            data[index_x] = math.sqrt(comp[0][index_x] ** 2 + comp[1][index_x] ** 2)

        return data

    def read_variable_barotropic_sea_water_from_direction_at_time(self, date):
        comp = self.read_variable_barotropic_sea_water_velocity_at_time(date)

        data = np.zeros([self.get_nb_points()])
        data[::] = np.NAN

        for index_x in range(0, self.get_nb_points()):
            data[index_x] = (180.0 + math.degrees(math.atan2(comp[1][index_x], comp[0][index_x]))) % 360.0

        return data

    def read_variable_barotropic_sea_water_to_direction_at_time(self, date):
        comp = self.read_variable_barotropic_sea_water_velocity_at_time(date)

        data = np.zeros([self.get_nb_points()])
        data[::] = np.NAN

        for index_x in range(0, self.get_nb_points()):
            data[index_x] = (math.degrees(math.atan2(comp[1][index_x], comp[0][index_x]))) % 360.0

        return data

    # METEO
    def read_variable_wind_10m_at_time(self, date):
        index_t = self.find_time_index(date)

        data = np.zeros([2, self.get_nb_points()])
        data[::] = np.NAN

        if len(index_t) > 1:
            layers = np.zeros([len(index_t),2, self.get_nb_points()])
            layers[::] = np.NAN

            for t in range(0, len(index_t)):
                comp = self.reader.read_variable_wind_10m_at_time(index_t[t])
                layers[t][0] = comp[0]
                layers[t][1] = comp[1]

            data[0] = self.interpolate_time(date,index_t,layers[:,0,:])
            data[1] = self.interpolate_time(date, index_t, layers[:,1,:])

        else:
            data = self.reader.read_variable_wind_10m_at_time(index_t[0])

        return data

    def read_variable_wind_speed_10m_at_time(self, date):
        comp = self.read_variable_wind_10m_at_time(date)

        data = np.zeros([self.get_nb_points()])
        data[::] = np.NAN

        for index_x in range(0,self.get_nb_points()):

            data[index_x] = math.sqrt(comp[0][index_x] ** 2 + comp[1][index_x] ** 2)

        return data

    def read_variable_wind_from_direction_10m_at_time(self, date):
        comp = self.read_variable_wind_10m_at_time(date)

        data = np.zeros([self.get_nb_points()])
        data[::] = np.NAN

        for index_x in range(0, self.get_nb_points()):
            data[index_x] = 270. - (math.degrees(math.atan2(comp[1][index_x], comp[0][index_x]))) + 180.0 % 360.0

        return data

    def read_variable_wind_to_direction_10m_at_time(self, date):
        comp = self.read_variable_wind_10m_at_time(date)

        data = np.zeros([self.get_nb_points()])
        data[::] = np.NAN

        for index_x in range(0, self.get_nb_points()):
            data[index_x] = 270. - (math.degrees(math.atan2(comp[1][index_x], comp[0][index_x]))) % 360.0

        return data

    def read_variable_surface_air_pressure_at_time(self, date):
        index_t = self.find_time_index(date)

        if len(index_t) > 1:
            layers = np.zeros([len(index_t), self.get_nb_points()])
            layers[::] = np.NAN

            for t in range(0, len(index_t)):
                layers[t] = self.reader.read_variable_surface_air_pressure_at_time(index_t[t])

            data = self.interpolate_time(date,index_t,layers)

        else:
            data = self.reader.read_variable_surface_air_pressure_at_time(index_t[0])

        return data

    def read_variable_rainfall_amount_at_time(self, date):
        index_t = self.find_time_index(date)

        if len(index_t) > 1:
            layers = np.zeros([len(index_t), self.get_nb_points()])
            layers[::] = np.NAN

            for t in range(0, len(index_t)):
                layers[t] = self.reader.read_variable_rainfall_amount_at_time(index_t[t])

            data = self.interpolate_time(date,index_t,layers)

        else:
            data = self.reader.read_variable_rainfall_amount_at_time(index_t[0])

        return data

    def read_variable_water_volume_transport_into_sea_water_from_rivers_at_time(self, date):
        index_t = self.find_time_index(date)

        if len(index_t) > 1:
            layers = np.zeros([len(index_t), self.get_nb_points()])
            layers[::] = np.NAN

            for t in range(0, len(index_t)):
                layers[t] = self.reader.read_variable_water_volume_transport_into_sea_water_from_rivers_at_time(index_t[t])

            data = self.interpolate_time(date, index_t, layers)

        else:
            data = self.reader.read_variable_water_volume_transport_into_sea_water_from_rivers_at_time(index_t[0])

        return data




    # TODO trouver une place propre pour ce qui suit

    def read_variable_sea_surface_temperature(self):

        data = self.reader.read_variable_sea_surface_temperature();
        return self.interpolate_all_times(data)

    def read_variable_sea_surface_salinity(self):

        data = self.reader.read_variable_sea_surface_salinity()
        return self.interpolate_all_times(data)

    def read_variable_sea_water_pressure_at_sea_water_surface(self):

        data = self.reader.read_variable_sea_water_pressure_at_sea_water_surface()
        return self.interpolate_all_times(data)

    def read_variable_sea_surface_density(self):

        data = self.reader.read_variable_sea_surface_density()
        return self.interpolate_all_times(data)

    def read_variable_sea_water_turbidity(self):

        data = self.reader.read_variable_sea_water_turbidity()
        return self.interpolate_all_times(data)

    def read_variable_sea_water_electrical_conductivity(self):

        data = self.reader.read_variable_sea_water_electrical_conductivity()
        return self.interpolate_all_times(data)




