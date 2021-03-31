# -*- encoding:utf-8 -*-
"""
==================================
NetCDF SYMPHONIE (:mod:`coverage.io.netcdf.symphonie`)
==================================

.. currentmodule:: coverage.io.netcdf.symphonie

"""
from __future__ import division, print_function, absolute_import


from .SYMPHONIEBathycoteInReader import SYMPHONIEBathycoteInReader
from .SYMPHONIEReader import SYMPHONIEReader
from spatialetl.coverage.io.netcdf.symphonie.v273.SymphonieOfflineReader import SymphonieOfflineReader

__all__ = ['SYMPHONIEBathycoteInReader', 'SYMPHONIEReader', 'SymphonieOfflineReader']
