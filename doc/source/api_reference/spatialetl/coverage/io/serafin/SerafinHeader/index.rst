spatialetl.coverage.io.serafin.SerafinHeader
============================================

.. py:module:: spatialetl.coverage.io.serafin.SerafinHeader


Attributes
----------

.. autoapisummary::

   spatialetl.coverage.io.serafin.SerafinHeader.SLF_EIT


Classes
-------

.. autoapisummary::

   spatialetl.coverage.io.serafin.SerafinHeader.SerafinHeader


Module Contents
---------------

.. py:data:: SLF_EIT
   :value: 'iso-8859-1'


.. py:class:: SerafinHeader(title='', format_type='SERAFIN ', lang='fr', endian='>')

   !
   @brief: Data type for reading and storing the Serafin file header

   # Attributes:
   ## Sizes
   - file_size <int>: total file size
   - header_size <int>: size taken by header only (set by `_set_header_size`)
   - frame_size <int>: size taken by a single frame (set by `_set_frame_size`)

   ## Float precision: simple or double precision
       (attributes are set by `_set_as_single_precision` or `_set_as_double_precision`)
   - float_type <str>: 'f' or 'd'
   - float_size <int>: 4 or 8
   - np_float_type: np.float32 or np.float64

   - endian <str>: file endianness ('>' = big-endian, '<' = little-endian)

   - title <bytes:72>: title of the simulation
   - file_format <bytes:8>: file format type (e.g. b'SERAFIN ')
   - language <str>: language for variables ('fr' or 'en')

   - is_2d <bool>: True in 2D else False
   - has_knolg <bool>: True if ipobo is replaced by KNOLG array
   - date <[int]>: list with 6 integers (year, month, day, hour, minute and second) or None if not available
   - nb_frames <int>: number of frames

   - nb_var <int>: number of variables
   - nb_var_quadratic <int>: ?
   - nb_nodes_per_elem <int>: number of nodes per element
   - nb_planes <int>: number of layers (should be 0 in 2D)

   - var_IDs <[str]>: variable identifiers
   - var_names <[str]>: variable names
   - var_units <[str]>: variable units

   - params <(int)>: tuple of 10 integer parameters
   - nb_elements <int>: number of 3D elements
   - nb_nodes <int>: number of 3D nodes
   - nb_nodes_per_elem <int>: number of nodes per element (= 3 in 2D and 6 in 3D)
   - nb_nodes_2d <int>: number of 2D nodes (equals to `nb_nodes` in 2D)

   - mesh_origin <(float, float)>: x and y shift to apply to written coordinates (set by `set_mesh_origin`)
   - x_stored <numpy.1D-array>: east written coordinates [shape = nb_nodes]
   - y_stored <numpy.1D-array>: north written coordinates [shape = nb_nodes]
   - x <numpy.1D-array>: east coordinates [shape = nb_nodes] (set by `_compute_mesh_coordinates`)
   - y <numpy.1D-array>: north coordinates [shape = nb_nodes] (set by `_compute_mesh_coordinates`)
   - ikle <numpy.1D-array>: connectivity table with 1-indexed nodes number [shape = nb_nodes_per_elem * nb_elements]
   - ikle_2d <numpy.2D-array>: reshaped connectivity table
       [shape = (nb_elements, nb_nodes_per_elem)] (set by `_build_ikle_2d`)
   - ipobo <numpy.1D-array>: array with 0 values for inner nodes and 1-indexed node number for boundary nodes
       [shape = nb_nodes]


   .. py:attribute:: file_size
      :value: -1



   .. py:attribute:: header_size
      :value: -1



   .. py:attribute:: frame_size
      :value: -1



   .. py:attribute:: float_type
      :value: ''



   .. py:attribute:: float_size
      :value: -1



   .. py:attribute:: np_float_type
      :value: None



   .. py:attribute:: endian
      :value: '>'



   .. py:attribute:: title
      :value: b''



   .. py:attribute:: file_format
      :value: b''



   .. py:attribute:: language
      :value: 'fr'



   .. py:attribute:: is_2d
      :value: True



   .. py:attribute:: has_knolg
      :value: False



   .. py:attribute:: date
      :value: None



   .. py:attribute:: nb_var
      :value: 0



   .. py:attribute:: nb_var_quadratic
      :value: 0



   .. py:attribute:: nb_planes
      :value: 0



   .. py:attribute:: nb_frames
      :value: 0



   .. py:attribute:: var_IDs
      :value: []



   .. py:attribute:: var_names
      :value: []



   .. py:attribute:: var_units
      :value: []



   .. py:attribute:: params
      :value: (1, 0, 0, 0, 0, 0, 0, 0, 0, 0)



   .. py:attribute:: nb_elements
      :value: -1



   .. py:attribute:: nb_nodes
      :value: -1



   .. py:attribute:: nb_nodes_per_elem
      :value: -1



   .. py:attribute:: nb_nodes_2d
      :value: -1



   .. py:attribute:: mesh_origin
      :value: (0, 0)



   .. py:attribute:: x_stored
      :value: None



   .. py:attribute:: y_stored
      :value: None



   .. py:attribute:: x
      :value: None



   .. py:attribute:: y
      :value: None



   .. py:attribute:: ikle
      :value: None



   .. py:attribute:: ikle_2d
      :value: None



   .. py:attribute:: ipobo
      :value: None



   .. py:method:: _check_dim()


   .. py:method:: _set_file_format_and_precision(file_format)

      !
      @brief: Set some attributes to read/write single or double precision Serafin depending on file format
      @file_format <str>: Serafin file format (max length is 8)
      If file format is not recognized, the file is expected to be simple precision



   .. py:method:: _set_has_knolg()


   .. py:method:: _build_ikle_2d()


   .. py:method:: _set_as_single_precision()

      Set Serafin as single precision



   .. py:method:: _set_as_double_precision()

      Set Serafin as double precision



   .. py:method:: unpack_int(bytes2unpack, nb=1)

      Unpack bytes to integer(s)
      @param bytes2unpack <bytes>: bytes to unpack
      @param nb <int>: number of consecutive integers
      @return <(int)>



   .. py:method:: unpack_float(bytes2unpack, nb=1)

      Unpack bytes to float(s)
      @param bytes2unpack <bytes>: bytes to unpack
      @param nb <int>: number of consecutive float
      @return <(float)>



   .. py:method:: pack_int(*args, nb=1)

      Pack integers
      @param args <(int)>: integers to pack
      @param nb <int>: number of integer values
      @return <bytes>



   .. py:method:: pack_float(*args, nb=1)

      Pack floats
      @param args <(float)>: floats to pack
      @param nb <int>: number of float values
      @return <bytes>



   .. py:method:: toggle_endianness()

      Toggle original endianness (between big or little endian)



   .. py:method:: _compute_mesh_coordinates()

      Compute mesh coordinates from origin



   .. py:method:: _set_header_size()

      Set header size



   .. py:method:: _set_frame_size()

      Set frame size (all variable values for one time step)



   .. py:method:: _expected_file_size()

      Returns expected file size



   .. py:method:: summary()


   .. py:method:: copy()

      Returns a deep copy of the current instance



   .. py:method:: copy_as_2d()

      !
      Returns a 2D equivalent copy of the current instance
      @return <slf.Serafin.SerafinHeader>: output Serafin header



   .. py:method:: copy_as_3d(nb_planes)

      !
      Returns a 3D equivalent copy of the current instance
      @param nb_planes <int>: number of planes
      @return <slf.Serafin.SerafinHeader>: output Serafin header



   .. py:method:: nearest_node(target_x, target_y)

      !
      Find the nearest node of a target point (from x and y coordinates)
      @param target_x <float>: east target coordinate
      @param target_y <float>: north target coordinate
      @return <int>: node number (1-indexed)



   .. py:method:: same_2d_mesh(other)

      !
      @brief: Check if the other mesh is strictly identical (same order of nodes and elements) on the horizontal
      @param other <SerafinHeader>: header to compare
      @return: bool



   .. py:method:: is_double_precision()


   .. py:method:: to_single_precision()


   .. py:method:: empty_variables()

      Empty all variables



   .. py:method:: add_variable(var_ID, var_name, var_unit)

      !
      @brief: Add a single variable
      @param var_ID <str>: variables identifier (abbreviation)
      @param var_name <bytes>: variable name
      @param var_unit <bytes>: variable unit



   .. py:method:: add_variable_str(var_ID, var_name, var_unit)

      !
      @brief: Add a single variable with strings
      @param var_ID <str>: variables identifier (abbreviation)
      @param var_name <str>: variable name
      @param var_unit <str>: variable unit



   .. py:method:: add_variable_from_ID(var_ID)

      !
      @brief: Add new variable from its ID
      @param var_ID <str>: variable ID



   .. py:method:: set_mesh_origin(x, y)

      !
      @brief: Set mesh origin coordinates
      @param x <float>: X-coordinate of the origin
      @param y <float>: Y-coordinate of the origin



   .. py:method:: set_variables(selected_vars)

      !
      @brief: Set new variables
      @param selected_vars <[str, bytes, bytes]>: list composed of variable ID, name and unit



   .. py:method:: iter_on_all_variables()

      !
      Iterate on variables
      @return <str, bytes, bytes>: variable ID, name and unit



   .. py:method:: transform_mesh_copy(transformations)

      !
      @brief Apply transformations on 2D nodes of a mesh copy
      @param transformations <[geom.transformation.Transformation]>: list of successive transformations
      @return: modified copy of original header with transformed mesh nodes



   .. py:method:: transform_mesh(transformations)

      !
      @brief Apply transformations on mesh nodes (only in 2D)
      @param transformations <[geom.transformation.Transformation]>: list of successive transformations



   .. py:method:: get_all_edges()

      Get all edges (pair of nodes)



   .. py:method:: get_external_edges()

      Get external edges (pair of nodes)



   .. py:method:: iter_on_boundaries()

      Iterate over all boundaries to get the nodes describing them
      Each element is an ordered list of nodes (1-indexed numbering)

      How boundaries are defined:
      - the first boundary is the external boundary and then the islands (internal boundaries) are described
          in an arbitrary order
      - the external boundary is described counter-clockwise
      - the islands are described clockwise
      - the first boundary node should be the node with the minimum value of x+y (corresponds to the bottom left corner)



   .. py:method:: build_ipobo()

      Build IPOBO array containing 0 values for inner nodes and 1-indexed node number for boundary nodes
      /!\ This method will probably crash if some nodes are duplicated!



   .. py:method:: _set_as_2d()

      !
      Beware: nb_nodes_2d has to be consistant



   .. py:method:: from_triangulation(nodes, ikle)

      !
      Set to Serafin 2D header from a given triangulation
      @param nodes <numpy 2D-array>: x and y coordinates
      @param ikle <numpy 2D-array>: connectivity table (1-indexed)



   .. py:method:: from_file(file, file_size)

      !
      @param file <_io.BufferedReader>: input Serafin stream
      @param file_size <int>: file size (in bytes)
      @return <slf.Serafin.SerafinHeader>: output Serafin header



