spatialetl.utils.bing.core.tilesystem
=====================================

.. py:module:: spatialetl.utils.bing.core.tilesystem

.. autoapi-nested-parse::

   __author__ = Linlin Chen
   __email__ = lchen96@hawk.iit.edu

   @Description:
   This module implements a set of static methods used for Bing maps tile system.

   reference: https://msdn.microsoft.com/en-us/library/bb259689.aspx



Classes
-------

.. autoapisummary::

   spatialetl.utils.bing.core.tilesystem.TileSystem


Module Contents
---------------

.. py:class:: TileSystem

   Bases: :py:obj:`object`


   .. py:attribute:: EARTHRADIUS
      :value: 6378137



   .. py:attribute:: MAXLEVEL
      :value: 23



   .. py:method:: clip(val, minval, maxval)
      :staticmethod:


      Clips a number to be specified minval and maxval values.

      :param val {[double]} -- [value to be clipped]:
      :param minval {[double]} -- [minimal value bound]:
      :param maxval {[double]} -- [maximal value bound]:

      :returns: [double] -- [clipped value]



   .. py:method:: map_size(level)
      :staticmethod:


      Determines the map width and height (in pixels) at a specified level

      :param level {[int]} -- [level of detail:
      :type level {[int]} -- [level of detail: lowest detail) to 23 (highest detail
      :param from 1:
      :type from 1: lowest detail) to 23 (highest detail

      :returns: [int] -- [The map width and height in pixels (width == height)]



   .. py:method:: ground_resolution(lat, level)
      :staticmethod:


      Determines the ground resolution (in meters per pixel) at a specified latitude and level of detail

      :param lat {[double]} -- [latitude in degrees at which to measure the ground resolution]:
      :param level {[int]} -- [level of detail:
      :type level {[int]} -- [level of detail: lowest detail) to 23 (highest detail
      :param from 1:
      :type from 1: lowest detail) to 23 (highest detail

      :returns: [double] -- [The ground resolution, in meters per pixel]



   .. py:method:: map_scale(lat, level, screenDpi)
      :staticmethod:


      Determines the map scale at a specified latitude, level of detail, and screen resolution

      :param lat {[double]} -- [latitude in degrees at which to measure the ground resolution]:
      :param level {[type]} -- [level of detail:
      :type level {[type]} -- [level of detail: lowest detail) to 23 (highest detail
      :param from 1:
      :type from 1: lowest detail) to 23 (highest detail
      :param screenDpi {[type]} -- [Resolution of the screen:
      :param in dots per inch]:

      :returns: N]
      :rtype: [double] -- [The map scale, expressed as the denominator N of the ratio 1



   .. py:method:: latlong_to_pixelXY(lat, long, level)
      :staticmethod:


      Converts a point from latitude/longitude WGS-84 coordinates (in degrees)
      into pixel XY coordinates a specified level of detail

      :param lat {[double]} -- [Latitude of the point:
      :param in degrees]:
      :param long {[double]} -- [Longitude of the point:
      :param in degrees]:
      :param level {[int]} -- [Level of detail:
      :type level {[int]} -- [Level of detail: lowest detail) to 23 (highest detail
      :param from 1:
      :type from 1: lowest detail) to 23 (highest detail

      :returns: [int, int] -- [X coordinates in pixels; Y coordinates in pixels]



   .. py:method:: pixelXY_to_latlong(pixelX, pixelY, level)
      :staticmethod:


      Converts a pixel from pixel XY coordinates at a specified level of detail
      into latitude/longitude WGS-84 coordinates (in degrees)

      :param pixelX {[int]} -- [X coordinate of the point:
      :param in pixels.]:
      :param pixelY {[int]} -- [Y coordinate of the point:
      :param in pixels.]:
      :param level {[int]} -- [Level of detail:
      :type level {[int]} -- [Level of detail: lowest detail) to 23 (highest detail
      :param from 1:
      :type from 1: lowest detail) to 23 (highest detail

      :returns: [double, double] -- [Latitude in degrees; Longitude in degrees]



   .. py:method:: pixelXY_to_tileXY(pixelX, pixelY)
      :staticmethod:


      Converts pixel XY coordinates into tile XY coordinates of the tile containing the specified pixel.

      :param pixelX {[int]} -- [Pixel X coordinate]:
      :param pixelY {[int]} -- [Pixel Y coordinate]:

      :returns: [int, int] -- [Tile X coordinate; Tile Y coordinate]



   .. py:method:: tileXY_to_pixelXY(tileX, tileY)
      :staticmethod:


      Converts tile XY coordinates into pixel XY coordinates of the upper-left pixel of the specified tile

      :param tileX {[int]} -- [Tile X coordinate]:
      :param tileY {[int]} -- [Tile Y coordinate]:

      :returns: [int, int] -- [pixel X coordinate; pixel Y coordinate]



   .. py:method:: tileXY_to_quadkey(tileX, tileY, level)
      :staticmethod:


      Converts tile XY coordinates into a QuadKey at a specified level of detail
      interleaving tileY with tileX

      :param tileX {[int]} -- [Tile X coordinate]:
      :param tileY {[int]} -- [Tile Y coordinate]:
      :param level {[int]} -- [Level of detail:
      :type level {[int]} -- [Level of detail: lowest detail) to 23 (highest detail
      :param from 1:
      :type from 1: lowest detail) to 23 (highest detail

      :returns: [string] -- [A string containing the QuadKey]



   .. py:method:: quadkey_to_tileXY(quadkey)
      :staticmethod:


      Converts a QuadKey into tile XY coordinate

      :param quadkey {[string]} -- [QuadKey of the tile]:

      :returns: [int, int] -- [Tile X coordinate; Tile Y coordinate]



