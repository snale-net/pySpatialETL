spatialetl.coverage.io.serafin.SerafinReader
============================================

.. py:module:: spatialetl.coverage.io.serafin.SerafinReader

.. autoapi-nested-parse::

   !
   Read/Write Serafin files and manipulate associated data.

   Handles Serafin file with:
     - single and double precision
     - big and little endian

   The sizes (file, header, frame) are in bytes (8 bits).



Classes
-------

.. autoapisummary::

   spatialetl.coverage.io.serafin.SerafinReader.SerafinReader


Module Contents
---------------

.. py:class:: SerafinReader(myFilename, language='fr')

   Bases: :py:obj:`spatialetl.coverage.io.CoverageReader.CoverageReader`


   .. py:attribute:: language
      :value: 'fr'



   .. py:attribute:: file


   .. py:attribute:: file_size


   .. py:attribute:: header


   .. py:attribute:: x_size


   .. py:attribute:: y_size


   .. py:attribute:: time_ref


   .. py:method:: close()


   .. py:method:: is_regular_grid()


   .. py:method:: get_x_size()


   .. py:method:: get_y_size()


   .. py:method:: get_t_size()


   .. py:method:: read_axis_x(xmin, xmax, ymin, ymax)


   .. py:method:: read_axis_y(xmin, xmax, ymin, ymax)


   .. py:method:: read_axis_t(tmin, tmax, timestamp)


   .. py:method:: read_variable_longitude(xmin, xmax, ymin, ymax)


   .. py:method:: read_variable_latitude(xmin, xmax, ymin, ymax)


   .. py:method:: read_variable_time(tmin, tmax)


   .. py:method:: _get_var_index(var_ID)

      !
      @brief Handle data request by variable ID
      @param var_ID <str>: the ID of the requested variable
      @return index <int> the index of the frame (0-based)



   .. py:method:: read_variable_bathymetry(xmin, xmax, ymin, ymax)


   .. py:method:: read_var_in_frame(time_index, var_ID)

      !
      @brief Read a single variable in a frame
      @param time_index <int>: the index of the frame (0-based)
      @param var_ID <str>: variable ID
      @return <numpy 1D-array>: values of the variables, of length equal to the number of nodes



   .. py:method:: read_var_in_frame_as_3d(time_index, var_ID)

      !
      @brief Read a single variable in a 3D frame
      @param time_index <int>: the index of the frame (0-based)
      @param var_ID <str>: variable ID
      @return <numpy 2D-array>: values of the variables with shape (planes number, number of 2D nodes)



   .. py:method:: read_var_in_frame_at_layer(time_index, var_ID, iplan)

      !
      @brief Read a single variable in a frame at specific layer
      @param time_index <int>: the index of the frame (0-based)
      @param var_ID <str>: variable ID
      @param iplan <int>: 1-based index of layer
      @return <numpy 1D-array>: values of the variables, of length equal to the number of nodes



