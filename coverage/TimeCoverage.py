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
from coverage.Coverage import Coverage
from datetime import timedelta
from datetime import datetime
import cftime
from numpy import int,int32,int64
import math
import numpy as np
from array_split import shape_split
from coverage.operator.interpolator.InterpolatorCore import resample_2d_to_grid
import logging

class TimeCoverage(Coverage):
    """La classe TimeCoverage est une extension de la classe Coverage.
Elle rajoute une dimension temporelle à la couverture horizontale classique.
    """
    
    TIME_DATUM = datetime(1970, 1, 1)    
    TIME_DELTA = timedelta(minutes = 15)
    TIME_OVERLAPING_SIZE = 0

    def __init__(self, myReader,bbox=None,resolution_x=None,resolution_y=None):
            
        Coverage.__init__(self,myReader,bbox=bbox,resolution_x=resolution_x,resolution_y=resolution_y);
        self.times = self.read_axis_t();
        self.global_source_t_size = len(self.times)

    def create_mpi_map(self):

        self.global_source_t_size = self.reader.get_t_size()

        self.map_mpi = np.empty(self.size, dtype=object)
        source_sample = np.zeros([self.global_source_t_size, self.global_source_y_size, self.global_source_x_size])
        target_sample = np.zeros([self.global_source_t_size, self.global_target_y_size, self.global_target_x_size])

        # Découpage sur 2 axes
        source_slices = shape_split(source_sample.shape, self.size, axis=[0, 0, 0])

        slice_index = 0
        for slyce in source_slices.flatten():
            slice = tuple(slyce)

            map = {}
            # Grille source
            map["source_global_t_min"] = slice[0].start
            map["source_global_t_max"] = slice[0].stop
            map["source_global_x_min"] = slice[2].start
            map["source_global_x_max"] = slice[2].stop
            map["source_global_y_min"] = slice[1].start
            map["source_global_y_max"] = slice[1].stop

            map["source_t_size"] = map["source_global_t_max"] - map["source_global_t_min"]
            map["source_x_size"] = map["source_global_x_max"] - map["source_global_x_min"]
            map["source_y_size"] = map["source_global_y_max"] - map["source_global_y_min"]

            map["source_global_t_min_overlap"] = max(0, map["source_global_t_min"] - TimeCoverage.TIME_OVERLAPING_SIZE)
            map["source_global_t_max_overlap"] = min(self.global_source_t_size, map["source_global_t_max"] + TimeCoverage.TIME_OVERLAPING_SIZE)
            map["source_global_x_min_overlap"] = max(0,map["source_global_x_min"] - Coverage.HORIZONTAL_OVERLAPING_SIZE)
            map["source_global_x_max_overlap"] = min(self.global_source_x_size,map["source_global_x_max"] + Coverage.HORIZONTAL_OVERLAPING_SIZE)
            map["source_global_y_min_overlap"] = max(0,map["source_global_y_min"] - Coverage.HORIZONTAL_OVERLAPING_SIZE)
            map["source_global_y_max_overlap"] = min(self.global_source_y_size,map["source_global_y_max"] + Coverage.HORIZONTAL_OVERLAPING_SIZE)

            map["source_t_size_overlap"] = map["source_global_t_max_overlap"] - map["source_global_t_min_overlap"]
            map["source_x_size_overlap"] = map["source_global_x_max_overlap"] - map["source_global_x_min_overlap"]
            map["source_y_size_overlap"] = map["source_global_y_max_overlap"] - map["source_global_y_min_overlap"]

            map["source_t_min"] = TimeCoverage.TIME_OVERLAPING_SIZE
            map["source_t_max"] = map["source_t_size_overlap"] - TimeCoverage.TIME_OVERLAPING_SIZE
            map["source_x_min"] = Coverage.HORIZONTAL_OVERLAPING_SIZE
            map["source_x_max"] = map["source_x_size_overlap"] - Coverage.HORIZONTAL_OVERLAPING_SIZE
            map["source_y_min"] = Coverage.HORIZONTAL_OVERLAPING_SIZE
            map["source_y_max"] = map["source_y_size_overlap"] - Coverage.HORIZONTAL_OVERLAPING_SIZE

            if map["source_global_t_min"] == 0:
                map["source_t_min"] = 0

            if map["source_global_t_max"] == self.global_source_t_size:
                map["source_t_max"] = map["source_t_size_overlap"]

            if map["source_global_x_min"] == 0:
                map["source_x_min"] = 0

            if map["source_global_y_min"] == 0:
                map["source_y_min"] = 0

            if map["source_global_x_max"] == self.global_source_x_size:
                map["source_x_max"] = map["source_x_size_overlap"]

            if map["source_global_y_max"] == self.global_source_y_size:
                map["source_y_max"] = map["source_y_size_overlap"]

            # Grille destination
            # lon = self.reader.read_axis_x(0, self.global_source_x_size, 0,
            #                               self.global_source_y_size)
            # lat = self.reader.read_axis_y(0, self.global_source_x_size, 0,
            #                               self.global_source_y_size)
            #
            #
            # lBcorner = self.find_point_index(lon[map["source_global_y_min"],map["source_global_x_min"]],lat[map["source_global_y_min"],map["source_global_x_min"]])
            # hRcorner = self.find_point_index(lon[map["source_global_y_max"]-1,map["source_global_x_max"]-1],lat[map["source_global_y_max"]-1,map["source_global_x_max"]-1])

            lon = self.reader.read_axis_x(map["source_global_x_min"], map["source_global_x_max"],
                                          map["source_global_y_min"],
                                          map["source_global_y_max"])
            lat = self.reader.read_axis_y(map["source_global_x_min"], map["source_global_x_max"],
                                          map["source_global_y_min"],
                                          map["source_global_y_max"])

            Ymin = np.min(lat)
            Ymax = np.max(lat)
            Xmin = np.min(lon)
            Xmax = np.max(lon)

            lBcorner = self.find_point_index(Xmin,Ymin)
            hRcorner = self.find_point_index(Xmax,Ymax)

            print(self.rank,"bottom left",lBcorner)
            print(self.rank,"up right",hRcorner)

            map["target_global_t_min"] = map["source_global_t_min"]
            map["target_global_t_max"] = map["source_global_t_max"]
            map["target_global_x_min"] = max(0,lBcorner[0]-1)
            map["target_global_x_max"] = min(self.global_target_x_size,hRcorner[0]+1)
            map["target_global_y_min"] = max(0,lBcorner[1]-1)
            map["target_global_y_max"] = max(self.global_target_y_size,hRcorner[1]+1)

            self.map_mpi[slice_index] = map

            slice_index = slice_index + 1
   
    # Axis
    def find_time_index(self,t):
        """Retourne l'index de la date la plus proche à TIME_DELTA_MIN prêt.
    @type t: datetime ou int
    @param t: date souhaitée ou l'index de la date souhaitée
    @return:  l'index de la date la plus proche à TIME_DELTA_MIN prêt ou une erreur si aucune date n'a pu être trouvée."""

        if type(t) == int or type(t) == int32 or type(t) == int64:

            if t < 0 or t >= self.get_t_size():
                raise ValueError("Time index have to range between 0 and " + str(
                    self.get_t_size() - 1) + ". Actually Time index = " + str(t))

            return t;

        if type(t) == datetime or type(t) == cftime._cftime.real_datetime:
            zero_delta = timedelta(minutes = 00)
            for i in range(0,self.get_t_size()):
                if t - self.times[i] == zero_delta or t - self.times[i] < TimeCoverage.TIME_DELTA:
                    return i

            raise ValueError(""+str(t)+" was not found. Maybe the TimeCoverage.TIME_DELTA_MIN ("+ str(TimeCoverage.TIME_DELTA)+") is too small or the date is out the range.")
        else:
            raise ValueError(""+str(t)+" have to be an integer or a datetime. Current type: "+str(type(t)))
    
    def read_axis_t(self,timestamp=0):
        """Retourne les valeurs de l'axe t.
    @param timestamp: égale 1 si le temps est souhaité en timestamp depuis TIME_DATUM.
    @return:  un tableau à une dimensions [z] au format datetime ou timestamp si timestamp=1."""
        return self.reader.read_axis_t(self.map_mpi[self.rank]["source_global_t_min_overlap"],
                                       self.map_mpi[self.rank]["source_global_t_max_overlap"],timestamp)[
            self.map_mpi[self.rank]["source_t_min"]:self.map_mpi[self.rank]["source_t_max"]];

    def get_global_t_size(self):
        return self.reader.get_t_size()
    
    def get_t_size(self):
        """Retourne la taille de l'axe t.
    @return:  un entier correspondant à la taille de l'axe t."""
        return self.global_source_t_size;
    
    # Variables
    def read_variable_2D_sea_binary_mask_at_time(self, t):
        """Retourne le masque à la date souhaitée sur toute la couverture horizontale.
    @type t: datetime ou l'index
    @param t: date souhaitée
    @return: un tableau en deux dimensions [y,x]."""

        index_t = self.find_time_index(t);

        return self.reader.read_variable_2D_sea_binary_mask_at_time(index_t)

    def read_variable_2D_wet_binary_mask_at_time(self, t):
        """Retourne le masque à la date souhaitée sur toute la couverture horizontale.
    @type t: datetime ou l'index
    @param t: date souhaitée
    @return: un tableau en deux dimensions [y,x]."""

        index_t = self.find_time_index(t);

        return self.reader.read_variable_2D_wet_binary_mask_at_time(index_t)

    #################
    # HYDRO
    # Sea Surface
    #################
    def read_variable_sea_surface_height_above_mean_sea_level_at_time(self,t):

        index_t = self.find_time_index(t);
       
        data = self.reader.read_variable_sea_surface_height_above_mean_sea_level_at_time(
            self.map_mpi[self.rank]["source_global_t_min_overlap"]+index_t,
            self.map_mpi[self.rank]["source_global_x_min_overlap"],
            self.map_mpi[self.rank]["source_global_x_max_overlap"],
            self.map_mpi[self.rank]["source_global_y_min_overlap"],
            self.map_mpi[self.rank]["source_global_y_max_overlap"])

        if self.horizontal_resampling:

            data = resample_2d_to_grid(self.read_axis_x(type="source", with_overlap=True),
                                       self.read_axis_y(type="source", with_overlap=True),
                                       self.read_axis_x(type="target", with_overlap=False),
                                       self.read_axis_y(type="target", with_overlap=False),
                                       data,
                                       Coverage.HORIZONTAL_INTERPOLATION_METHOD)

            return data

        else:
            return data[
                   self.map_mpi[self.rank]["source_y_min"]:self.map_mpi[self.rank]["source_y_max"],
                   self.map_mpi[self.rank]["source_x_min"]:self.map_mpi[self.rank]["source_x_max"]]

    def read_variable_sea_surface_height_above_geoid_at_time(self, t):

        index_t = self.find_time_index(t);

        return self.reader.read_variable_sea_surface_height_above_geoid_at_time(index_t)

    def read_variable_sea_surface_temperature_at_time(self, t):
        """Retourne la temperature de surface à la date souhaitée
    @type t: datetime ou l'index
    @param t: date souhaitée
    @return: un tableau en deux dimensions [y,x]."""

        index_t = self.find_time_index(t);

        return self.reader.read_variable_sea_surface_temperature_at_time(index_t)

    def read_variable_sea_surface_salinity_at_time(self, t):
        """Retourne la salinité de surface à la date souhaitée
    @type t: datetime ou l'index
    @param t: date souhaitée
    @return: un tableau en deux dimensions [y,x]."""

        index_t = self.find_time_index(t);

        return self.reader.read_variable_sea_surface_salinity_at_time(index_t)

    def read_variable_sea_surface_pressure_at_time(self, t):
        """Retourne la pression à la surface de la mer (sea surface pressure) à la date souhaitée sur toute la couverture horizontale.
    @type t: datetime ou l'index
    @param t: date souhaitée
    @return: un tableau en deux dimensions [y,x]."""

        index_t = self.find_time_index(t);

        return self.reader.read_variable_sea_surface_pressure_at_time(index_t)

    def read_variable_sea_surface_density_at_time(self, t):
        """Retourne la densité de l'eau de surface à la date souhaitée
    @type t: datetime ou l'index
    @param t: date souhaitée
    @return: un tableau en deux dimensions [y,x]."""

        index_t = self.find_time_index(t);

        return self.reader.read_variable_sea_surface_density_at_time(index_t)

    def read_variable_sea_water_turbidity_at_time(self, t):
        """Retourne la turbidité de l'eau de surface à la date souhaitée
    @type t: datetime ou l'index
    @param t: date souhaitée
    @return: un tableau en deux dimensions [y,x]."""

        index_t = self.find_time_index(t);

        return self.reader.read_variable_sea_water_turbidity_at_time(index_t)


    def read_variable_sea_water_velocity_at_sea_water_surface_at_time(self, t):
        """Retourne les composantes u,v du courant à la date souhaitée
    @type t: datetime ou l'index
    @param t: date souhaitée
    @return: un tableau en deux dimensions [u_comp,v_comp] contenant chacun deux dimensions [y,x]."""

        index_t = self.find_time_index(t);

        return self.reader.read_variable_sea_water_velocity_at_sea_water_surface_at_time(index_t)

    #################
    # HYDRO
    # Ground level
    #################

    def read_variable_sea_water_temperature_at_ground_level_at_time(self, t):
        """Retourne la temperature de surface à la date souhaitée
    @type t: datetime ou l'index
    @param t: date souhaitée
    @return: un tableau en deux dimensions [y,x]."""

        index_t = self.find_time_index(t);

        return self.reader.read_variable_sea_water_temperature_at_ground_level_at_time(index_t)

    def read_variable_sea_water_salinity_at_ground_level_at_time(self, t):
        """Retourne la salinité de surface à la date souhaitée
    @type t: datetime ou l'index
    @param t: date souhaitée
    @return: un tableau en deux dimensions [y,x]."""

        index_t = self.find_time_index(t);

        return self.reader.read_variable_sea_water_salinity_at_ground_level_at_time(index_t)

    def read_variable_sea_water_velocity_at_ground_level_at_time(self, t):
        """Retourne les composantes u,v du courant à la date souhaitée
           @type t: datetime ou l'index
           @param t: date souhaitée
           @return: un tableau en deux dimensions [u_comp,v_comp] contenant chacun deux dimensions [y,x]."""

        index_t = self.find_time_index(t);

        return self.reader.read_variable_sea_water_velocity_at_ground_level_at_time(index_t)

    #################
    # HYDRO
    # 2D
    #################
    def read_variable_barotropic_sea_water_velocity_at_time(self,t):
        """Retourne les composantes u,v du courant à la date souhaitée
    @type t: datetime ou l'index
    @param t: date souhaitée
    @return: un tableau en deux dimensions [u_comp,v_comp] contenant chacun deux dimensions [y,x]."""

        index_t = self.find_time_index(t);

        return self.reader.read_variable_barotropic_sea_water_velocity_at_time(index_t)

    #################
    # WAVES
    # Sea Surface
    #################
    def read_variable_sea_surface_wave_significant_height_at_time(self,t):
        """Retourne la hauteur significative des vagues à la date souhaitée sur toute la couverture horizontale.
    @type t: datetime ou l'index
    @param t: date souhaitée
    @return: un tableau en deux dimensions [y,x]."""

        index_t = self.find_time_index(t);

        return self.reader.read_variable_sea_surface_wave_significant_height_at_time(index_t)

    def read_variable_sea_surface_wave_breaking_height_at_time(self, t):
        """Retourne la hauteur de déferlement des vagues à la date souhaitée sur toute la couverture horizontale.
    @type t: datetime ou l'index
    @param t: date souhaitée
    @return: un tableau en deux dimensions [y,x]."""

        index_t = self.find_time_index(t);

        return self.reader.read_variable_sea_surface_wave_breaking_height_at_time(index_t)
    
    def read_variable_sea_surface_wave_mean_period_at_time(self,t):

        index_t = self.find_time_index(t);

        return self.reader.read_variable_sea_surface_wave_mean_period_at_time(index_t)

    def read_variable_sea_surface_wave_peak_period_at_time(self,t):

        index_t = self.find_time_index(t);

        return self.reader.read_variable_sea_surface_wave_peak_period_at_time(index_t)

    def read_variable_sea_surface_wave_from_direction_at_time(self,t):

        index_t = self.find_time_index(t);

        return self.reader.read_variable_sea_surface_wave_from_direction_at_time(index_t)

    def read_variable_sea_surface_wave_to_direction_at_time(self,t):

        index_t = self.find_time_index(t);

        return self.reader.read_variable_sea_surface_wave_to_direction_at_time(index_t)

    def read_variable_sea_surface_wave_stokes_drift_velocity_at_time(self, t):
        """Retourne la dérive de Stokes en surface à la date souhaitée sur toute la couverture horizontale.
    @type t: datetime ou l'index
    @param t: date souhaitée
    @return: un tableau en deux dimensions [y,x]."""

        index_t = self.find_time_index(t);

        return self.reader.read_variable_sea_surface_wave_stokes_drift_velocity_at_time(index_t)

    def read_variable_radiation_pressure_bernouilli_head_at_time(self, t):
        """Retourne la pression J due aux vagues à la date souhaitée sur toute la couverture horizontale.
    @type t: datetime ou l'index
    @param t: date souhaitée
    @return: un tableau en deux dimensions [y,x]."""

        index_t = self.find_time_index(t);

        return self.reader.read_variable_radiation_pressure_bernouilli_head_at_time(index_t)

    def read_variable_sea_surface_wave_energy_flux_to_ocean_at_time(self, t):
        """Retourne la waves_to_ocean_energy_flux à la date souhaitée sur toute la couverture horizontale.
    @type t: datetime ou l'index
    @param t: date souhaitée
    @return: un tableau en deux dimensions [y,x]."""

        index_t = self.find_time_index(t);

        return self.reader.read_variable_sea_surface_wave_energy_flux_to_ocean_at_time(index_t)

    def read_variable_sea_surface_wave_energy_dissipation_at_ground_level_at_time(self, t):
        """Retourne la l'énergie des vagues dissipée par le fond à la date souhaitée sur toute la couverture horizontale.
    @type t: datetime ou l'index
    @param t: date souhaitée
    @return: un tableau en deux dimensions [y,x]."""

        index_t = self.find_time_index(t);

        return self.reader.read_variable_sea_surface_wave_energy_dissipation_at_ground_level_at_time(index_t)


    #################
    # WAVES
    # Momentum flux
    #################
    def read_variable_atmosphere_momentum_flux_to_waves_at_time(self,t):
        """Retourne la composante u du tau atmosphere->vagues à la date souhaitée sur toute la couverture horizontale.
    @type t: datetime ou l'index
    @param t: date souhaitée
    @return: un tableau en deux dimensions [y,x]."""

        index_t = self.find_time_index(t);

        return self.reader.read_variable_atmosphere_momentum_flux_to_waves_at_time(index_t)
    
    def read_variable_waves_momentum_flux_to_ocean_at_time(self,t):
        """Retourne la composante u du tau vagues->ocean à la date souhaitée sur toute la couverture horizontale.
    @type t: datetime ou l'index
    @param t: date souhaitée
    @return: un tableau en deux dimensions [y,x]."""

        index_t = self.find_time_index(t);

        return self.reader.read_variable_waves_momentum_flux_to_ocean_at_time(index_t)

    #################
    # METEO
    # 2D
    #################
    def read_variable_rainfall_amount_at_time(self, t):
        """Retourne les composantes u,v de rain à la date souhaitée
    @type t: datetime ou l'index
    @param t: date souhaitée
    @return: un tableau en deux dimensions [u_comp,v_comp] contenant chacun deux dimensions [y,x]."""

        index_t = self.find_time_index(t);

        return self.reader.read_variable_rainfall_amount_at_time(index_t)

    #################
    # METEO
    # Sea surface
    #################

    def read_variable_surface_air_pressure_at_time(self,t):
        """Retourne la pression à la surface à la date souhaitée sur toute la couverture horizontale.
    @type t: datetime ou l'index
    @param t: date souhaitée
    @return: un tableau en deux dimensions [y,x]."""

        index_t = self.find_time_index(t);

        return self.reader.read_variable_surface_air_pressure_at_time(index_t)

    def read_variable_sea_surface_air_pressure_at_time(self, t):
        """Retourne la pression à la surface à la date souhaitée sur toute la couverture horizontale.
    @type t: datetime ou l'index
    @param t: date souhaitée
    @return: un tableau en deux dimensions [y,x]."""

        index_t = self.find_time_index(t);

        return self.reader.read_variable_sea_surface_air_pressure_at_time(index_t)

    def read_variable_wind_stress_at_time(self, t):
        """Retourne les composantes u,v de la contrainte de vent à la date souhaitée
    @type t: datetime ou l'index
    @param t: date souhaitée
    @return: un tableau en deux dimensions [u_comp,v_comp] contenant chacun deux dimensions [y,x]."""

        index_t = self.find_time_index(t);

        return self.reader.read_variable_wind_stress_at_time(index_t)

    def read_variable_surface_downward_sensible_heat_flux_at_time(self, t):
        """Retourne les composantes u,v de surface sensible heat flux à la date souhaitée
    @type t: datetime ou l'index
    @param t: date souhaitée
    @return: un tableau en deux dimensions [u_comp,v_comp] contenant chacun deux dimensions [y,x]."""

        index_t = self.find_time_index(t);

        return self.reader.read_variable_surface_downward_sensible_heat_flux_at_time(index_t)

    def read_variable_surface_downward_latent_heat_flux_at_time(self, t):
        """Retourne les composantes u,v de surface latente heat flux à la date souhaitée
    @type t: datetime ou l'index
    @param t: date souhaitée
    @return: un tableau en deux dimensions [u_comp,v_comp] contenant chacun deux dimensions [y,x]."""

        index_t = self.find_time_index(t);

        return self.reader.read_variable_surface_downward_latent_heat_flux_at_time(index_t)

    def read_variable_surface_air_temperature_at_time(self, t):
        """Retourne les composantes u,v de surface air temperature à la date souhaitée
    @type t: datetime ou l'index
    @param t: date souhaitée
    @return: un tableau en deux dimensions [u_comp,v_comp] contenant chacun deux dimensions [y,x]."""

        index_t = self.find_time_index(t);

        return self.reader.read_variable_surface_air_temperature_at_time(index_t)

    def read_variable_dew_point_temperature_at_time(self, t):
        """Retourne les composantes u,v de dewpoint temperature à la date souhaitée
    @type t: datetime ou l'index
    @param t: date souhaitée
    @return: un tableau en deux dimensions [u_comp,v_comp] contenant chacun deux dimensions [y,x]."""

        index_t = self.find_time_index(t);

        return self.reader.read_variable_dew_point_temperature_at_time(index_t)

    def read_variable_surface_downwards_solar_radiation_at_time(self, t):
        """Retourne les composantes u,v de surface solar radiation downwards à la date souhaitée
    @type t: datetime ou l'index
    @param t: date souhaitée
    @return: un tableau en deux dimensions [u_comp,v_comp] contenant chacun deux dimensions [y,x]."""

        index_t = self.find_time_index(t);

        return self.reader.read_variable_surface_downwards_solar_radiation_at_time(index_t)

    def read_variable_surface_downwards_thermal_radiatio_at_time(self, t):
        """Retourne les composantes u,v de surface thermal radiation downwards à la date souhaitée
    @type t: datetime ou l'index
    @param t: date souhaitée
    @return: un tableau en deux dimensions [u_comp,v_comp] contenant chacun deux dimensions [y,x]."""

        index_t = self.find_time_index(t);

        return self.reader.read_variable_surface_downwards_thermal_radiatio_at_time(index_t)

    def read_variable_surface_solar_radiation_at_time(self, t):
        """Retourne les composantes u,v de surface solar radiation à la date souhaitée
    @type t: datetime ou l'index
    @param t: date souhaitée
    @return: un tableau en deux dimensions [u_comp,v_comp] contenant chacun deux dimensions [y,x]."""

        index_t = self.find_time_index(t);

        return self.reader.read_variable_surface_solar_radiation_at_time(index_t)

    def read_variable_surface_thermal_radiation_at_time(self, t):
        """Retourne les composantes u,v de surface thermal radiation à la date souhaitée
    @type t: datetime ou l'index
    @param t: date souhaitée
    @return: un tableau en deux dimensions [u_comp,v_comp] contenant chacun deux dimensions [y,x]."""

        index_t = self.find_time_index(t);

        return self.reader.read_variable_surface_thermal_radiation_at_time(index_t)

    #################
    # METEO
    # At 10 m
    #################
    def read_variable_wind_10m_at_time(self,t):
        """Retourne les composantes u,v du vent à la date souhaitée
    @type t: datetime ou l'index
    @param t: date souhaitée
    @return: un tableau en deux dimensions [u_comp,v_comp] contenant chacun deux dimensions [y,x]."""

        index_t = self.find_time_index(t);

        return self.reader.read_variable_wind_10m_at_time(index_t)

    def read_variable_wind_speed_10m_at_time(self, date):
        comp = self.read_variable_wind_10m_at_time(date)
        result = np.zeros([self.get_y_size(),self.get_x_size()])
        result[:] = np.nan
        for x in range(0,self.get_x_size()):
            for y in range(0,self.get_y_size()):
                result[y,x] = math.sqrt(comp[0][y,x] ** 2 + comp[1][y,x] ** 2)

        return result

    def read_variable_wind_from_direction_10m_at_time(self, date):
        comp = self.read_variable_wind_10m_at_time(date)
        result = np.zeros([self.get_y_size(), self.get_x_size()])
        result[:] = np.nan
        for x in range(0, self.get_x_size()):
            for y in range(0, self.get_y_size()):
                result[y, x] = 270. - (180.0 / math.pi) * (math.atan2(comp[0][y,x], comp[1][y,x])) + 180.0 % 360.0

        return result

    def read_variable_wind_to_direction_10m_at_time(self, date):
        comp = self.read_variable_wind_10m_at_time(date)
        result = np.zeros([self.get_y_size(), self.get_x_size()])
        result[:] = np.nan
        for x in range(0, self.get_x_size()):
            for y in range(0, self.get_y_size()):
                result[y, x] = 270. - (180.0 / math.pi) * (math.atan2(comp[0][y,x], comp[1][y,x])) % 360.0

        return result









