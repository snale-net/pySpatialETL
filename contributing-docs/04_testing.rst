Testing Infrastructure
======================

.. contents:: On this page
   :local:
   :depth: 2

Unit Tests
----------

Location: Tests follow the code structure
~~~~~~~~~

.. code-block:: text

   spatialetl-core/spatialetl/coverage/
   ├── coverage.py
   └── tests/
       └── test_coverage.py

Execution
~~~~~~~~~

.. code-block:: bash

   # Test a specific module
   uv run pytest spatialetl-core/spatialetl/coverage/tests/

   # Run all tests
   uv run pytest

**Requirement**: Every PR must include at least one unit test, except for documentation-only changes.

Writing Tests
-------------

Example Test
~~~~~~~~~~~~

.. code-block:: python

   # spatialetl-core/spatialetl/coverage/tests/test_coverage.py
   import numpy as np
   import pytest
   from spatialetl.coverage.coverage import Coverage
   from spatialetl.coverage.io.memory_reader import MemoryReader


   def test_bathymetry_simple_grid():
       """Test bathymetry reading on a 5x5 grid."""
       # Prepare test data
       test_data = np.zeros([5, 5])
       np.fill_diagonal(test_data, 2)

       reader = MemoryReader(x=range(5), y=range(5), bathy=test_data)
       coverage = Coverage(reader=reader, nb_thread=2)

       # Verify the result
       result = coverage.read_variable_bathymetry(thread_id=0)
       np.testing.assert_array_equal(result, test_data)


   def test_invalid_thread_count():
       """Test thread count validation."""
       reader = MemoryReader(x=range(5), y=range(5), bathy=np.zeros([5, 5]))

       with pytest.raises(ValueError):
           Coverage(reader=reader, nb_thread=0)

Best Practices
~~~~~~~~~~~~~~

- **Clear naming**: ``test_<function>_<scenario>``
- **One test = one behavior**: Keep each test focused on a single aspect
- **Test edge cases**: null values, expected errors, empty grids
- **Use pytest fixtures** to share setup code across tests