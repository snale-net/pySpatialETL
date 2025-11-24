Contributors' Guide
===================

Welcome to the pySpatialETL Contributors' Guide!

pySpatialETL is an ETL library for geospatial data processing. This guide will help you get started with contributing to the project.

Documentation Structure
-----------------------

This guide is organized into the following sections:

- `Quick Start Guide <01_quick_start.rst>`_ -
- `Local Development Setup <02_local_virtualenv.rst>`_ - Detailed environment setup with UV
- `Pull Request Guidelines <03_pull_requests.rst>`_ - How to submit PRs
- `Testing Guide <04_testing.rst>`_ - Running and writing tests

Project Overview
----------------

**Project Structure**

.. code-block:: text

   pySpatialETL/
   ├── spatialetl-core/          # Core library
   │   └── spatialetl/
   │       ├── coverage/         # Gridded data processing
   │       ├── exception/
   │       ├── point/            # Multi-point data
   │       ├── operator/         # Grid operators
   │       └── utils/            # Utilities
   │
   ├── providers/                # Data format providers
   │   ├── common/
   │   │   ├── gdal/
   │   │   ├── grib/
   │   │   ├── mpi/
   │   │   └── netcdf/
   │   ├── ecmwf/
   │   ├── gmt/
   │   ├── hycom/
   │   ├── mercator/
   │   ├── meteofrance/
   │   ├── swan/
   │   ├── symphonie/
   │   ├── telemac/
   │   └── ww3
   ├── contributings-docs/
   ├── demo/                     # Example scripts
   └── doc/                      # Documentation

Need Help?
----------

1. Check the documentation in this directory
2. Search existing issues on GitHub
3. Open a new issue or discussion