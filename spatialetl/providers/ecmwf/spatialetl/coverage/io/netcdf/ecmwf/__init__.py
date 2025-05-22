# -*- encoding:utf-8 -*-
"""
==================================
NetCDF ECMWF (:mod:`coverage.io.netcdf.ecmwf`)
==================================

.. currentmodule:: coverage.io.netcdf.ecmwf

"""
from __future__ import division, print_function, absolute_import

from .ecmwf_reader import ECMWFReader
from .ecmwf_writer import ECMWFWriter

__all__ = ['ECMWFReader', 'ECMWFWriter']
