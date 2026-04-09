Providers
=========

pySpatialETL uses a modular provider architecture where each data source
is handled by an independent package.

.. toctree::
   :maxdepth: 1
   :caption: Common Utilities

   /api_reference/spatialetl/providers/common/gdal/index
   /api_reference/spatialetl/providers/common/grib/index
   /api_reference/spatialetl/providers/common/mpi/index
   /api_reference/spatialetl/providers/common/netcdf/index

.. toctree::
   :maxdepth: 1
   :caption: Data Source Providers

   /api_reference/spatialetl/providers/ecmwf/index
   /api_reference/spatialetl/providers/gmt/index
   /api_reference/spatialetl/providers/hycom/index
   /api_reference/spatialetl/providers/mercator/index
   /api_reference/spatialetl/providers/meteofrance/index
   /api_reference/spatialetl/providers/swan/index
   /api_reference/spatialetl/providers/symphonie/index
   /api_reference/spatialetl/providers/telemac/index
   /api_reference/spatialetl/providers/ww3/index

Installation
------------

Each provider is installed as an optional dependency using UV:

.. code-block:: bash

   # Install specific provider
   uv sync --package spatialetl-providers-meteofrance

   # Install all providers with core
   uv sync