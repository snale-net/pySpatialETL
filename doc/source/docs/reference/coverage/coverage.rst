Class Coverage
==============
Lorem ipsum dolor sit amet, consectetur adipiscing elit. Sed do eiusmod tempor incididunt ut labore et dolore magna aliqua. Ut enim ad minim veniam, quis nostrud exercitation ullamco laboris nisi ut aliquip ex ea commodo consequat. Duis aute irure dolor in reprehenderit in voluptate velit esse cillum dolore eu fugiat nulla pariatur. Excepteur sint occaecat cupidatat non proident, sunt in culpa qui officia deserunt mollit anim id est laborum.

.. py:module:: spatialetl.coverage.coverage

Classes
-------

.. autoapisummary::

   spatialetl.coverage.coverage.Coverage

Static attributes
-----------------

.. autosummary::
   :toctree: _generated
   :nosignatures:

   Coverage.HORIZONTAL_INTERPOLATION_METHOD
   Coverage.HORIZONTAL_OVERLAPING_SIZE

Class definition
----------------

.. autoclass:: Coverage
   :noindex:
   :show-inheritance:
   :exclude-members: HORIZONTAL_INTERPOLATION_METHOD, HORIZONTAL_OVERLAPING_SIZE, __init__

Constructeur
------------

.. automethod:: Coverage.__init__
   :noindex:

Methods
-------

.. autosummary::
   :toctree: _generated
   :nosignatures:

   Coverage.check_bbox_validity
   Coverage.check_point_is_inside
   Coverage.create_mpi_map
   Coverage.update_mpi_map
   Coverage.read_metadata
   Coverage.get_x_size
   Coverage.get_y_size
   Coverage.is_regular_grid
   Coverage.read_axis_x
   Coverage.read_axis_y
   Coverage.find_point_index
   Coverage.read_variable_bathymetry
   Coverage.read_variable_topography
   Coverage.read_variable_mesh_size
   Coverage.read_variable_x_mesh_size
   Coverage.read_variable_y_mesh_size
   Coverage.read_variable_2D_sea_binary_mask
   Coverage.read_variable_Ha