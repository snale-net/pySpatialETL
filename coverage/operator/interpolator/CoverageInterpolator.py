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

from coverage.operator.interpolator.InterpolatorCore import resample_2d_to_grid
import numpy as np

class CoverageInterpolator():
    """
    Cette classe permet d'interpoler un coverage non régulier sur une grille régulière.
    """

    def __init__(self, cov,resX,resY,depths):
        """
    Constructeur
    @param cov : la coverage
    @param resX : résolution souhaitée en X, en degrès
    @param resY : résolution souhaitée en Y, en degrès
    """
        self.coverage = cov
        self.targetResX = resX
        self.targetResY = resY;

        # we compute the destination grid
        Ymin=np.min(self.coverage.read_axis_y())
        Ymax=np.max(self.coverage.read_axis_y())
        Xmin=np.min(self.coverage.read_axis_x())
        Xmax=np.max(self.coverage.read_axis_x())

        res=np.mean([self.targetResX,self.targetResY])
        self.lon_reg,self.lat_reg=np.meshgrid(np.arange(Xmin, Xmax, res),np.arange(Ymin, Ymax, res))

        self.depths = depths

    # Axis
    def read_axis_x(self):
        return self.lon_reg[0, :]

    def read_axis_y(self):
        return self.lat_reg[:,0]

    def read_axis_z(self):
        return self.depths

    def read_axis_t(self,timestamp=0):
        return self.coverage.read_axis_t(timestamp)

    #Scalar
    # HYDRO
    def read_variable_2D_sea_binary_mask(self):
        return resample_2d_to_grid(self.coverage.read_axis_x(),self.coverage.read_axis_y(),self.lon_reg,self.lat_reg,self.coverage.read_variable_2D_sea_binary_mask())

    def read_variable_2D_wet_binary_mask_at_time(self, time):
        return resample_2d_to_grid(self.coverage.read_axis_x(), self.coverage.read_axis_y(), self.lon_reg, self.lat_reg,
                                   self.coverage.read_variable_2D_wet_binary_mask_at_time(time))

    def read_variable_topography(self):
        return resample_2d_to_grid(self.coverage.read_axis_x(),self.coverage.read_axis_y(),self.lon_reg,self.lat_reg,self.coverage.read_variable_topography())

    def read_variable_bathymetry(self):

        return resample_2d_to_grid(self.coverage.read_axis_x(),self.coverage.read_axis_y(),self.lon_reg,self.lat_reg,self.coverage.read_variable_bathymetry())

    def read_variable_mesh_size(self):
        return resample_2d_to_grid(self.coverage.read_axis_x(),self.coverage.read_axis_y(),self.lon_reg,self.lat_reg,self.coverage.read_variable_mesh_size())

    def read_variable_Ha(self):
        return resample_2d_to_grid(self.coverage.read_axis_x(),self.coverage.read_axis_y(),self.lon_reg,self.lat_reg,self.coverage.read_variable_Ha())


    def read_variable_sea_surface_height_above_mean_sea_level_at_time(self,time):
        return resample_2d_to_grid(self.coverage.read_axis_x(),self.coverage.read_axis_y(),self.lon_reg,self.lat_reg,self.coverage.read_variable_sea_surface_height_above_mean_sea_level_at_time(time))

    def read_variable_sea_surface_height_above_geoid_at_time(self, time):
        return resample_2d_to_grid(self.coverage.read_axis_x(), self.coverage.read_axis_y(), self.lon_reg, self.lat_reg,
                                   self.coverage.read_variable_sea_surface_height_above_geoidl_at_time(time))

    def read_variable_bathy_ssh_at_time(self,time):
       return resample_2d_to_grid(self.coverage.read_axis_x(),self.coverage.read_axis_y(),self.lon_reg,self.lat_reg,self.coverage.read_variable_bathy_ssh_at_time(time))

    def read_variable_baroclinic_sea_water_velocity_at_time_and_depth(self,time,depth):

        cur = self.coverage.read_variable_baroclinic_sea_water_velocity_at_time_and_depth(time,depth)
        ucur = resample_2d_to_grid(self.coverage.read_axis_x(),self.coverage.read_axis_y(),self.lon_reg,self.lat_reg,cur[0])
        vcur = resample_2d_to_grid(self.coverage.read_axis_x(),self.coverage.read_axis_y(),self.lon_reg,self.lat_reg,cur[1])

        return [ucur,vcur]

    def read_variable_barotropic_sea_water_velocity_at_time(self,time):

        cur = self.coverage.read_variable_barotropic_sea_water_velocity_at_time(time)
        ucur = resample_2d_to_grid(self.coverage.read_axis_x(),self.coverage.read_axis_y(),self.lon_reg,self.lat_reg,cur[0])
        vcur = resample_2d_to_grid(self.coverage.read_axis_x(),self.coverage.read_axis_y(),self.lon_reg,self.lat_reg,cur[1])

        return [ucur,vcur]

    def read_variable_sea_water_salinity_at_time_and_depth(self, time, depth):
        return resample_2d_to_grid(self.coverage.read_axis_x(), self.coverage.read_axis_y(), self.lon_reg, self.lat_reg,
                                   self.coverage.read_variable_sea_water_salinity_at_time_and_depth(time, depth))

    def read_variable_sea_water_temperature_at_time_and_depth(self, time, depth):
        return resample_2d_to_grid(
            self.coverage.read_axis_x(), self.coverage.read_axis_y(), self.lon_reg, self.lat_reg,
            self.coverage.read_variable_sea_water_temperature_at_time_and_depth(time, depth))

    # WAVES
    def read_variable_sea_surface_wave_significant_height_at_time(self, time):
        return resample_2d_to_grid(self.coverage.read_axis_x(), self.coverage.read_axis_y(), self.lon_reg, self.lat_reg,
                                   self.coverage.read_variable_sea_surface_wave_significant_height_at_time(time))

    def read_variable_sea_surface_wave_breaking_height_at_time(self, time):
        return resample_2d_to_grid(self.coverage.read_axis_x(), self.coverage.read_axis_y(), self.lon_reg, self.lat_reg,
                                   self.coverage.read_variable_sea_surface_wave_breaking_height_at_time(time))

    def read_variable_sea_surface_wave_mean_period_at_time(self, time):
        return resample_2d_to_grid(self.coverage.read_axis_x(), self.coverage.read_axis_y(), self.lon_reg, self.lat_reg,
                                   self.coverage.read_variable_sea_surface_wave_mean_period_at_time(time))

    def read_variable_sea_surface_wave_peak_period_at_time(self, time):
        return resample_2d_to_grid(self.coverage.read_axis_x(), self.coverage.read_axis_y(), self.lon_reg, self.lat_reg,
                                   self.coverage.read_variable_sea_surface_wave_peak_period_at_time(time))

    def read_variable_sea_surface_wave_to_direction_at_time(self, time):
        return resample_2d_to_grid(self.coverage.read_axis_x(), self.coverage.read_axis_y(), self.lon_reg, self.lat_reg,
                                   self.coverage.read_variable_sea_surface_wave_to_direction_at_time(time))

    def read_variable_sea_surface_wave_from_direction_at_time(self, time):
        return resample_2d_to_grid(self.coverage.read_axis_x(), self.coverage.read_axis_y(), self.lon_reg, self.lat_reg,
                                   self.coverage.read_variable_sea_surface_wave_from_direction_at_time(time))

    def read_variable_sea_surface_wave_energy_dissipation_at_ground_level_at_time(self, time):
        return resample_2d_to_grid(self.coverage.read_axis_x(), self.coverage.read_axis_y(), self.lon_reg, self.lat_reg,
                                   self.coverage.read_variable_sea_surface_wave_energy_dissipation_at_ground_level_at_time(time))

    def read_variable_sea_surface_wave_energy_flux_to_ocean_at_time(self, time):
        return resample_2d_to_grid(self.coverage.read_axis_x(), self.coverage.read_axis_y(), self.lon_reg, self.lat_reg,
                                   self.coverage.read_variable_sea_surface_wave_energy_flux_to_ocean_at_time(
                                       time))

    def read_variable_atmosphere_momentum_flux_to_waves_at_time(self,time):

        cur = self.coverage.read_variable_atmosphere_momentum_flux_to_waves_at_time(time)
        ucur = resample_2d_to_grid(self.coverage.read_axis_x(),self.coverage.read_axis_y(),self.lon_reg,self.lat_reg,cur[0])
        vcur = resample_2d_to_grid(self.coverage.read_axis_x(),self.coverage.read_axis_y(),self.lon_reg,self.lat_reg,cur[1])

        return [ucur,vcur]

    def read_variable_waves_momentum_flux_to_ocean_at_time(self,time):

        cur = self.coverage.read_variable_waves_momentum_flux_to_ocean_at_time(time)
        ucur = resample_2d_to_grid(self.coverage.read_axis_x(),self.coverage.read_axis_y(),self.lon_reg,self.lat_reg,cur[0])
        vcur = resample_2d_to_grid(self.coverage.read_axis_x(),self.coverage.read_axis_y(),self.lon_reg,self.lat_reg,cur[1])

        return [ucur,vcur]

    def read_variable_sea_surface_wave_stokes_drift_velocity_at_time(self, time):
        cur = self.coverage.read_variable_sea_surface_wave_stokes_drift_velocity_at_time(time)
        ucur = resample_2d_to_grid(self.coverage.read_axis_x(), self.coverage.read_axis_y(), self.lon_reg, self.lat_reg,
                                   cur[0])
        vcur = resample_2d_to_grid(self.coverage.read_axis_x(), self.coverage.read_axis_y(), self.lon_reg, self.lat_reg,
                                   cur[1])
        return [ucur, vcur]

    def read_variable_radiation_pressure_bernouilli_head_at_time(self, time):
        return resample_2d_to_grid(self.coverage.read_axis_x(), self.coverage.read_axis_y(), self.lon_reg, self.lat_reg,
                                   self.coverage.read_variable_radiation_pressure_bernouilli_head_at_time(time))


    # METEO
    def read_variable_wind_10m_at_time(self,time):

        cur = self.coverage.read_variable_wind_10m_at_time(time)
        ucur = resample_2d_to_grid(self.coverage.read_axis_x(),self.coverage.read_axis_y(),self.lon_reg,self.lat_reg,cur[0])
        vcur = resample_2d_to_grid(self.coverage.read_axis_x(),self.coverage.read_axis_y(),self.lon_reg,self.lat_reg,cur[1])

        return [ucur,vcur]

    def read_variable_wind_stress_at_time(self,time):

        cur = self.coverage.read_variable_wind_stress_at_time(time)
        ucur = resample_2d_to_grid(self.coverage.read_axis_x(),self.coverage.read_axis_y(),self.lon_reg,self.lat_reg,cur[0])
        vcur = resample_2d_to_grid(self.coverage.read_axis_x(),self.coverage.read_axis_y(),self.lon_reg,self.lat_reg,cur[1])
        return [ucur,vcur]