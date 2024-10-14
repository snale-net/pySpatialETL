#! /usr/bin/env python2.7
# -*- coding: utf-8 -*-
import io
import os
from setuptools import setup, find_packages

def read(*paths, **kwargs):
    """Read the contents of a text file safely.
    >>> read("project_name", "VERSION")
    '0.1.0'
    >>> read("README.md")
    ...
    """

    content = ""
    with io.open(
            os.path.join(os.path.dirname(__file__), *paths),
            encoding=kwargs.get("encoding", "utf8"),
    ) as open_file:
        content = open_file.read().strip()
    return content


def read_requirements(path):
    return [
        line.strip()
        for line in read(path).split("\n")
        if not line.startswith(('"', "#", "-", "git+"))
    ]


setup (
       name='pySpatialETL',
       version=read("VERSION"),
       packages=find_packages(),

       # Declare your packages' dependencies here, for eg:
       install_requires=read_requirements("requirements.txt"),

       author='Fabien Retif',
       author_email='fabien.retif@snale.net',

       #summary = 'Just another Python package for the cheese shop',
       url='https://lab.snale.net/produits/pyspatialetl',
       license='',
       long_description=' Extract - Transform - Load (ETL) of geospatial data, that is, digital geospatial information representing space and time-varying phenomena. This python package deals with multi-dimensional gridded data and multi-dimensional multi-point data.',

       # could also include long_description, download_url, classifiers, etc.
       )