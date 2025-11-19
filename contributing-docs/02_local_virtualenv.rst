Local Development with UV
==========================

This guide describes how to install pySpatialETL in development mode with UV.

Why UV?
-------

UV is a fast, modern alternative to pip and virtualenv that:

- **Faster** - 10-100x faster than traditional tools
- **Reliable** - Uses lockfiles for reproducible installations
- **Modern** - Built for Python 3.9+ with modern packaging standards
- **Simple** - Unified tool for environments, packages, and dependencies

For more information, consult the `official UV documentation <https://docs.astral.sh/uv/>`_.

Prerequisites & Clone the Repository
-------------

Ensure you have Git and UV installed.  see the `Quick Start Guide <01_quick_start.rst>`_ for installation instructions.

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

Complete installation (all providers):

.. code-block:: bash

   uv sync --all-extras

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

   uv run pytest

Working on Core Only
--------------------

When you only want to work on spatialetl-core, you can run ``uv sync`` in the ``spatialetl-core`` folder:

.. code-block:: bash

   cd spatialetl-core
   uv sync

This will install all dependencies needed to run tests for the core library.

Working on Individual Providers
--------------------------------

Sometimes you want to only work on a specific provider and you only want to install that provider's dependencies and run only that provider's tests.

This can be done very easily with ``uv`` by going to the provider's folder and running ``uv sync`` there. For example, to install dependencies of the ``netcdf`` provider:

.. code-block:: bash

   cd providers/netcdf
   uv sync

This will use the ``.venv`` environment in the root of your project and will install dependency of your provider and providers it depends on and its development dependencies.

Then running tests for the provider is as simple as:

.. code-block:: bash

   uv run pytest

Note that the ``uv sync`` command will automatically synchronize all dependencies needed for your provider and its development dependencies.

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
   uv run pytest

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

Corrupted Environment
~~~~~~~~~~~~~~~~~~~~~

In case of major issues, recreate the environment:

.. code-block:: bash

   rm -rf .venv
   uv venv
   uv sync

Quick Reference Commands
------------------------

.. code-block:: bash

   # Create environment
   uv venv

   # Install core only
   uv sync

   # Install with providers
   uv sync --extra netcdf --extra ecmwf

   # Install all providers
   uv sync --all-extras

   # Work on core only
   cd spatialetl-core && uv sync

   # Work on specific provider
   cd providers/netcdf && uv sync

   # List packages
   uv pip list

   # Run tests
   uv run pytest

   # Run a script
   uv run python script.py

   # Clean and restart
   rm -rf .venv && uv venv && uv sync
