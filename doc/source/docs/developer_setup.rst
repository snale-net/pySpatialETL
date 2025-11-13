Developer Setup
===============

This guide describes how to install pySpatialETL in development mode with UV.

Prerequisites
-------------

Before starting, ensure you have:

- UV installed (For more information, consult the `official UV documentation <https://docs.astral.sh/uv/>`_.)
- Git configured
- Python 3.9 or higher

Clone the Repository
--------------------

.. code-block:: bash

   git clone https://github.com/snale-net/pySpatialETL.git
   cd pySpatialETL
   git checkout develop

Installation with UV
--------------------

Create the Environment
~~~~~~~~~~~~~~~~~~~~~~

.. code-block:: bash

   uv venv

This command creates a virtual environment in ``.venv/``.

Install Core and Providers
~~~~~~~~~~~~~~~~~~~~~~~~~~~

Basic installation (core only):

.. code-block:: bash

   uv sync

Installation with specific providers:

.. code-block:: bash

   uv sync --extra netcdf --extra ecmwf --extra hycom

Complete installation (all available providers):

.. code-block:: bash

   uv sync --extra netcdf --extra ecmwf --extra hycom --extra mercator --extra meteofrance --extra swan --extra symphonie --extra telemac --extra ww3

Install pytest
~~~~~~~~~~~~~~

To run tests:

.. code-block:: bash

   uv pip install pytest

Verify Installation
-------------------

List Installed Packages
~~~~~~~~~~~~~~~~~~~~~~~~

.. code-block:: bash

   uv pip list | grep spatialetl

You should see ``spatialetl-core`` and the installed providers with their local paths.

Run Tests
~~~~~~~~~

.. code-block:: bash

   uv run pytest spatialetl-core/spatialetl/coverage/tests/


Daily Development
-----------------

Modify Code
~~~~~~~~~~~

Changes in ``spatialetl-core/`` or ``providers/`` are immediately reflected (editable mode).

Run Commands
~~~~~~~~~~~~

Always use ``uv run`` to execute commands in the environment:

.. code-block:: bash

   uv run python my_script.py
   uv run pytest tests/

Add a Dependency
~~~~~~~~~~~~~~~~

1. Modify the ``pyproject.toml`` file of the relevant package
2. Synchronize the environment:

.. code-block:: bash

   uv sync

Troubleshooting
---------------

"ModuleNotFoundError" Error
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

**Symptom**: ``ModuleNotFoundError: No module named 'XXX'``

**Solution**: Verify the package is installed and synchronize:

.. code-block:: bash

   uv sync
   uv pip list

Provider Build Error
~~~~~~~~~~~~~~~~~~~~

**Symptom**: ``Failed to build spatialetl-providers-XXX``

**Known cause**: The ``common-mpi`` and ``gmt`` providers have configuration issues (November 2025).

**Solution**: Install without these providers:

.. code-block:: bash

   uv sync --extra netcdf --extra ecmwf --extra hycom --extra mercator --extra meteofrance --extra swan --extra symphonie --extra telemac --extra ww3

Corrupted Environment
~~~~~~~~~~~~~~~~~~~~~

In case of major issues, recreate the environment:

.. code-block:: bash

   rm -rf .venv
   uv venv
   uv sync

Pytest Cannot Find Modules
~~~~~~~~~~~~~~~~~~~~~~~~~~~

Install pytest in the UV environment:

.. code-block:: bash

   uv pip install pytest

Known Limitations
-----------------

Non-functional Providers
~~~~~~~~~~~~~~~~~~~~~~~~

As of November 2025, the following providers have configuration errors:

- ``spatialetl-providers-common-mpi``: incorrect flit_core configuration
- ``spatialetl-providers-gmt``: empty [tool.flit.module] section

These providers will be fixed in a future release.


Quick Reference Commands
------------------------

.. code-block:: bash

   # Create environment
   uv venv

   # Install core only
   uv sync

   # Install with providers
   uv sync --extra netcdf --extra ecmwf

   # List packages
   uv pip list

   # Run tests
   uv run pytest

   # Run a script
   uv run python script.py

   # Clean and restart
   rm -rf .venv && uv venv && uv sync