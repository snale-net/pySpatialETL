#! /usr/bin/env python2.7
# -*- coding: utf-8 -*-

from setuptools import setup, find_packages

setup (
       name='pyGeoSpatialETL',
       version='0.1',
       packages=find_packages(),

       # Declare your packages' dependencies here, for eg:
       install_requires=['netCDF4>=1.0.7','scipy>=0.16.1','pandas>=0.17.1'],

       author='Fabien Retif',
       author_email='fabien.retif@zoho.com',

       #summary = 'Just another Python package for the cheese shop',
       url='',
       license='',
       long_description=' Extract - Transform - Load (ETL) of geospatial data, that is, digital geospatial information representing space and time-varying phenomena. This python package deals with multi-dimensional gridded data and multi-dimensional multi-point data.',

       # could also include long_description, download_url, classifiers, etc.
  
       )