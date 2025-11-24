Quick Start Guide
=================

This guide will help you make your first contribution to pySpatialETL in just a few steps.

Prerequisites
-------------

- Git installed and configured
- UV installed (see `UV documentation <https://docs.astral.sh/uv/>`_)

If you need help with UV installation:

.. code-block:: bash

   # macOS/Linux
   curl -LsSf https://astral.sh/uv/install.sh | sh

   # Windows (PowerShell)
   powershell -c "irm https://astral.sh/uv/install.ps1 | iex"

Step 1: Fork and Clone
-----------------------

1. Fork the repository on GitHub: https://github.com/snale-net/pySpatialETL

2. Clone your fork:

   .. code-block:: bash

      git clone https://github.com/YOUR_USERNAME/pySpatialETL.git
      cd pySpatialETL

3. Checkout the develop branch:

   .. code-block:: bash

      git checkout develop

Step 2: Set Up Environment
---------------------------

1. Create a virtual environment:

   .. code-block:: bash

      uv venv

   Or with a specific Python version:

   .. code-block:: bash

      uv venv --python 3.10.7

2. Install the package in development mode:

   Basic installation (core only):

   .. code-block:: bash

      uv sync

   Or with specific providers:

   .. code-block:: bash

      uv sync --extra netcdf --extra ecmwf --extra hycom

   Or install all providers:

   .. code-block:: bash

      uv sync --all-extras

Step 3: Verify Installation
----------------------------

Run tests to check everything works:

   .. code-block:: bash

     uv run pytest spatialetl-core/spatialetl/coverage/tests/

   All tests:

   .. code-block:: bash
     uv run pytest

If tests pass, you're ready to contribute!

Step 4: Make Your Changes
--------------------------

1. Create a new branch:

   .. code-block:: bash

      git checkout -b feature/my-feature

   Use descriptive branch names:

   - ``feature/add-gfs-reader`` for new features
   - ``fix/netcdf-time-handling`` for bug fixes
   - ``docs/improve-readme`` for documentation

2. Make your changes in ``spatialetl-core/`` or ``providers/``

3. Write or update tests

4. Run tests:

   .. code-block:: bash

      uv run pytest

Step 5: Submit Pull Request
----------------------------

1. Stage and commit your changes:

   .. code-block:: bash

      git add .
      git commit -m "Add support for GFS GRIB2 format"

2. Push to your fork:

   .. code-block:: bash

      git push -u origin feature/my-feature

3. Create a Pull Request on GitHub:

   - Go to https://github.com/snale-net/pySpatialETL
   - Click "Compare & pull request"
   - **Target branch:** ``develop`` (not ``master``!)
   - Fill in the PR description

Quick Reference Commands
------------------------

.. code-block:: bash

   # Clone and setup
   git clone https://github.com/YOUR_USERNAME/pySpatialETL.git
   cd pySpatialETL
   git checkout develop
   uv venv
   uv sync --extra netcdf

   # Create branch and make changes
   git checkout -b feature/my-feature
   # ... make your changes ...
   uv run pytest

   # Commit and push
   git add .
   git commit -m "Description of changes"
   git push -u origin feature/my-feature

Next Steps
----------

- Detailed setup guide: `Local Development Setup <02_local_virtualenv.rst>`_
- PR best practices: `Pull Request Guidelines <03_pull_requests.rst>`_
- Testing guide: `Testing Guide <04_testing.rst>`_