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
import numpy as np
import math
 
def distance_on_unit_sphere(long1, lat1,long2, lat2):
    """
    Calcule la distance en kilomètre entre les deux point. Les coordonnées
    sont données en longitude, latitude.

    @type  long1: number
    @param long1: Coordonnée X du point 1.
    @type  lat1: number
    @param lat1: Coordonnée Y du point 1.
    @type  long2: number
    @param long2: Coordonnée X du point 1.
    @type  lat2: number
    @param lat2: Coordonnée Y du point 1.
    @return:  la ditance en kilomètre en ces deux point sur la Terre.
    """

    # Convert latitude and longitude to 
    # spherical coordinates in radians.
    degrees_to_radians = math.pi/180.0
         
    # phi = 90 - latitude
    phi1 = (90.0 - lat1)*degrees_to_radians
    phi2 = (90.0 - lat2)*degrees_to_radians
         
    # theta = longitude
    theta1 = long1*degrees_to_radians
    theta2 = long2*degrees_to_radians
         
    # Compute spherical distance from spherical coordinates.
         
    # For two locations in spherical coordinates 
    # (1, theta, phi) and (1, theta', phi')
    # cosine( arc length ) = 
    #    sin phi sin phi' cos(theta-theta') + cos phi cos phi'
    # distance = rho * arc length
     
    cos = (math.sin(phi1)*math.sin(phi2)*math.cos(theta1 - theta2) + 
           math.cos(phi1)*math.cos(phi2))
    arc = math.acos( cos )
 
    # Remember to multiply arc by the radius of the earth 
    # in your favorite set of units to get length.
    return arc*6373

class Coverage(object):
    """
La classe Coverage représente une couverture spatiale sur l'horizontale. Les point qui représentent cette couverture
peuvent être alignés sur une maille régulière (x,y) ou sur une maille non-régulière ((x1,y1),(x2,y2)). En fonction du
type de maille, les fonctions de lecture des axes retourneront des tableaux à une ou deux dimensions. Pour éviter
un chargement en mémoire de la totalité du fichier, la coverage contient un pointeur vers un lecteur. Les couches
sont donc lues à la demande dans le fichier.

Attention, les axes sont toujours inversés dans les tableaux à cause de NetCDF.
Soit l'axe y en premier puis l'axe x. Exemple : [y,x]

@param  myReader: lecteur de fichier
"""
    def __init__(self, myReader):          
            
        self.reader = myReader;
        self.regular_grid = None

        # try to fill metadata
        self.read_metadata()

    # Read metadata
    def read_metadata(self):
        """
        Lit la metadonnée du fichier si le lecteur contient une fonction read_metadata()
        Returns
        -------

        """

        if  "read_metadata" in dir(self.reader):
            m = self.reader.read_metadata()

    # Axis        
    def read_axis_x(self):
        """Retourne les valeurs (souvent la longitude) de l'axe x.
    @return:  un tableau à une ou deux dimensions selon le type de maille des valeurs de l'axe x (souvent la longitude) : [x] ou [y,x]."""

        return self.reader.read_axis_x()
        
    def read_axis_y(self):
        """Retourne les valeurs (souvent la latitude) de l'axe y.
    @return:  un tableau à une ou deux dimensions selon le type de maille des valeurs de l'axe y (souvent la latitude) : [x] ou [y,x]."""
        return self.reader.read_axis_y()
    
    def get_x_size(self):
        """Retourne la taille de l'axe x.
    @return:  un entier correspondant à la taille de l'axe x."""
        if self.is_regular_grid():
            return np.shape(self.read_axis_x())[0];
        else:
            return np.shape(self.read_axis_x())[1];
    
    def get_y_size(self):
        """Retourne la taille de l'axe y.
    @return:  un entier correspondant à la taille de l'axe y."""
        if self.is_regular_grid():
            return np.shape(self.read_axis_y())[0];
        else:
            return np.shape(self.read_axis_y())[0];
    
    def is_regular_grid(self):
        """Retourne vrai si la maille est régulière, sinon faux.
    @return:  vrai si la maille est régulière sinon faux."""

        if self.regular_grid is None:
            if self.read_axis_x().ndim == 1:
                self.regular_grid = True
            else:
                self.regular_grid = False

        return self.regular_grid
        
    def find_point_index(self,target_lon,target_lat,method="classic"):
        """Retourne le point le plus proche du point donné en paramètre.
    @param target_lon: Coordonnée longitude du point
    @param target_lat: Coordonnée latitude du point
    @param method : Méthode de calcul. "Classic" = On parcourt toute la grille à la recherche du plus prêt.
    @return: un tableau contenant
     [0] : l'index x du point le plus proche
     [1] : l'index y du point le plus proche
     [2] : la coordonnée en longitude du point le plus proche
     [3] : la coordonnée en latitude point le plus proche
     [4] : la distance du point le plus proche en kilomètre."""
        min_dist=10000000
        lon = self.read_axis_x()
        lat = self.read_axis_y()
        mask = self.read_variable_2D_sea_binary_mask()
        
        nearest_x_index = 0
        nearest_y_index = 0
        nearest_lon = np.nan
        nearest_lat = np.nan

        if method=="classic":
            for x in range(0, self.get_x_size()):
                for y in range(0, self.get_y_size()):

                    if(mask[y,x] == 1): #=Terre

                        if self.is_regular_grid():
                            dist = distance_on_unit_sphere(target_lon,target_lat,lon[x],lat[y])
                        else:
                            dist = distance_on_unit_sphere(target_lon,target_lat,lon[y,x],lat[y,x])


                        if dist < min_dist:
                            min_dist = dist
                            if self.is_regular_grid():
                                nearest_lon = lon[x]
                                nearest_lat = lat[y]
                            else:
                                nearest_lon = lon[y,x]
                                nearest_lat = lat[y,x]

                            nearest_y_index = y
                            nearest_x_index = x
        else:
            raise RuntimeError("Method "+str(method)+" is not implemented yet.")

        if nearest_lon == np.nan or nearest_lat == np.nan:
            raise RuntimeError("No nearest point found.")
		
        return [nearest_x_index,nearest_y_index,nearest_lon,nearest_lat,min_dist]
    
    # Variables
    #################
    # HYDRO
    # 2D
    #################
    def read_variable_bathymetry(self):     
        """Retourne la bathymétrie sur toute la couverture
    @return: un tableau en deux dimensions [y,x]."""
        return self.reader.read_variable_bathymetry()
    
    def read_variable_mesh_size(self):     
        """Retourne la taille de la grille sur toute la couverture
    @return: un tableau en deux dimensions [y,x]."""
        return self.reader.read_variable_mesh_size()
    
    def read_variable_2D_sea_binary_mask(self):
        """Retourne le masque terre/mer sur toute la couverture
    @return: un tableau en deux dimensions [y,x].
            0 = Terre
            1 = Mer
    """
        return self.reader.read_variable_2D_sea_binary_mask()

    #################
    # METEO
    # 2D
    #################
    def read_variable_topography(self):     
        """Retourne la topographie sur toute la couverture
    @return: un tableau en deux dimensions [y,x]."""
        return self.reader.read_variable_topography()

    #OTHERS
    def read_variable_Ha(self):
        """Retourne l'amplitude de la réanalyse
    @return: un tableau en deux dimensions [y,x]."""
        return self.reader.read_variable_Ha()

    # def read_variable_sea_surface_density(self):
    #     """Retourne la densité de l'eau de surface
    #         @return: un tableau en deux dimensions [y,x]."""
    #     return self.reader.read_variable_sea_surface_density()
    #
    # def read_variable_sea_water_turbidity(self):
    #     """Retourne la turbidité de l'eau de surface
    #                 @return: un tableau en deux dimensions [y,x]."""
    #     return self.reader.read_variable_sea_water_turbidity()

        
        