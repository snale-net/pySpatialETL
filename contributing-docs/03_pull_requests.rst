Pull Request Guidelines
========================

This document explains how to create Pull Requests and details the code standards expected during their implementation.

.. contents:: On this page
   :local:
   :depth: 2

Before Submitting
-----------------

Pre-submission Checklist
~~~~~~~~~~~~~~~~~~~~~~~~

Before submitting a PR from your fork, ensure it meets the following requirements:

✔️ **Tests Required**

Include tests (doctests, unit tests with pytest, or both).
Tests must cover new features and edge cases.

✔️ **ReadTheDocs Build Passing**

Documentation build must pass without errors.
Maintainers will never merge a PR that breaks linting or documentation.

✔️ **All Conversations Resolved**

All discussions must be resolved before the PR can be merged.

✔️ **Rebase Over Merge**

Rebase your PR frequently to maintain a clean history and facilitate review.
All conflicts must be resolved.

✔️ **Squash and Merge**

Regardless of commit count during review, the PR will be merged as a single commit.
Maintainers may request you to clean up or consolidate commits before merging.

✔️ **MIT License Required**

All new files must begin with the MIT license header.

✔️ **Code + Tests + Docs in Same PR**

When adding a feature, documentation must be updated in the same PR:

- Docstrings in the code
- Sphinx documentation if applicable
- README files if relevant

✔️ **Small and Focused PRs**

Do not mix refactoring with new features.
Small PRs are easier to review and facilitate cherry-picking for patch releases.

For large changes:

1. Create a draft PR for discussion
2. Submit multiple smaller derived PRs

✔️ **Run Tests Locally**

Tests follow the same structure as the code.

Example: changes in ``spatialetl-core/spatialetl/coverage/`` ⇒ tests in ``spatialetl-core/spatialetl/coverage/tests/``

Run tests with UV:

.. code-block:: bash

   # Test a specific module
   uv run pytest spatialetl-core/spatialetl/coverage/tests/

   # Run all tests
   uv run pytest

✔️ **Test on Python 3.9**

Minimum supported version: Python 3.9.
Some recent features (match/case, new type hints) are not available in this version.

✔️ **Conventional Commit Messages**

Recommended format: ``[Type] Short description``

Accepted types:

- ``[Feat]``: New feature
- ``[Fix]``: Bug fix
- ``[Docs]``: Documentation only
- ``[Chore]``: Maintenance, dependencies
- ``[Refactor]``: Refactoring without functional changes
- ``[Test]``: Adding or modifying tests

Example: ``[Feat] add AROME forecast data provider``

Review Process
--------------

Conversation Resolution
~~~~~~~~~~~~~~~~~~~~~~~

A PR is mergeable only when **all conversations are resolved**.

This ensures:

- Clear visibility of PR status
- Faster review/merge cycle
- Limited use of "Request changes" to truly blocking issues

What Reviewers Expect
~~~~~~~~~~~~~~~~~~~~~

- Readable and well-documented code
- Relevant passing tests
- Up-to-date documentation

After Review
~~~~~~~~~~~~

If changes are requested:

1. Apply the corrections
2. Commit and push to your branch
3. Reply to comments indicating completion

Coding Style and Best Practices
--------------------------------

Don't Use Asserts Outside Tests
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Our community agreed that for various reasons we do not use ``assert`` in production
code of SpatialETL. Assertions are disabled when Python runs in optimized mode
(``python -O``), making them unreliable for production validation.

**❌ Avoid:**

.. code-block:: python

   assert grid_shape == (1800, 1536)
   assert geom.is_valid()

**✅ Prefer:**

.. code-block:: python

   if grid_shape != (1800, 1536):
       raise ValueError(f"Invalid AROME grid shape: expected (1800, 1536), got {grid_shape}")

   if not geom.is_valid():
       raise InvalidGeometryError(f"Geometry validation failed: {geom.wkt}")

**Exception:** Type checking guards are acceptable:

.. code-block:: python

   if TYPE_CHECKING:
       assert isinstance(dataset, xr.Dataset)

Use Standard Python Exceptions
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

We prioritize Python's standard exceptions over custom ones
for better interoperability and clarity.

**Hierarchy:**

1. Python standard exceptions (``ValueError``, ``TypeError``, ``KeyError``, ``OSError``)
2. Custom exceptions in ``spatialetl/exceptions.py`` (only when semantically necessary)

**❌ Too generic:**

.. code-block:: python

   raise SpatialETLException("Invalid coordinates")

**✅ Specific and standard:**

.. code-block:: python

   raise ValueError(f"Latitude must be between -90 and 90, got {lat}")
   raise KeyError(f"Missing required config key: {key}")

**✅ Custom when justified:**

.. code-block:: python

   # In spatialetl/exceptions.py
   class InvalidCRSError(ValueError):
       """Raised when coordinate reference system is invalid or unsupported."""
       pass

   # Usage
   if not crs.is_valid:
       raise InvalidCRSError(f"Unsupported CRS: {crs}")

Don't Use time() for Duration Calculations
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Use monotonic clocks for duration calculations to avoid issues with system time
adjustments (NTP synchronization, daylight saving time).

If you wish to compute the time difference between two events within the same process,
use ``time.monotonic()``, not ``time.time()`` nor ``datetime.now()``.

**❌ Unreliable (affected by NTP, DST):**

.. code-block:: python

   start = time.time()
   process_meteorological_data()
   duration = time.time() - start

**✅ Monotonic clock:**

.. code-block:: python

   start = time.monotonic()
   process_meteorological_data()
   duration = time.monotonic() - start
   log.info(f"Processing took {duration:.2f}s")

** High-precision benchmarking:**

If you are measuring duration for performance reasons, use ``time.perf_counter()``.
On many platforms, this uses the same underlying clock mechanism as monotonic, but
``perf_counter()`` is guaranteed to be the highest accuracy clock on the system.

.. code-block:: python

   start = time.perf_counter()
   extract_arome_grid()
   duration = time.perf_counter() - start

**⚠️ Database timestamps:** If the start time of a duration calculation needs to be
stored in a database, then this has to be done using ``datetime`` objects. In all
other cases, using ``datetime`` for duration calculation MUST be avoided as creating
and diffing datetime operations are (comparatively) slow.

Documentation Standards
~~~~~~~~~~~~~~~~~~~~~~~

**Docstrings**
   English, following `Google style <https://google.github.io/styleguide/pyguide.html#38-comments-and-docstrings>`_.

**Type hints**
   Required for all public APIs.

**Example:**

.. code-block:: python

   def extract_arome_variable(
       dataset: xr.Dataset,
       variable: str,
       bounds: tuple[float, float, float, float],
   ) -> np.ndarray:
       """Extract a meteorological variable from AROME dataset.

       Args:
           dataset: NetCDF dataset containing AROME forecast data.
           variable: Variable name (e.g., 'temperature', 'precipitation').
           bounds: Spatial bounds as (min_lon, min_lat, max_lon, max_lat).

       Returns:
           Extracted data array with spatial subset.

       Raises:
           KeyError: If variable not found in dataset.
           ValueError: If bounds are invalid.
       """
       if not is_valid_france_bounds(bounds):
           raise ValueError(f"Bounds outside France extent: {bounds}")

       return dataset[variable].sel(
           lon=slice(bounds[0], bounds[2]),
           lat=slice(bounds[1], bounds[3])
       ).values