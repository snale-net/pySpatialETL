spatialetl.coverage.io.ascii.swan.SWANUnstructuredReader
========================================================

.. py:module:: spatialetl.coverage.io.ascii.swan.SWANUnstructuredReader


Classes
-------

.. autoapisummary::

   spatialetl.coverage.io.ascii.swan.SWANUnstructuredReader.SWANUnstructuredReader
   spatialetl.coverage.io.ascii.swan.SWANUnstructuredReader.Converters


Functions
---------

.. autoapisummary::

   spatialetl.coverage.io.ascii.swan.SWANUnstructuredReader.swantime2datetime
   spatialetl.coverage.io.ascii.swan.SWANUnstructuredReader.dir2cat
   spatialetl.coverage.io.ascii.swan.SWANUnstructuredReader.deg2uv


Module Contents
---------------

.. py:class:: SWANUnstructuredReader(myNodefile, myFilename)

   Bases: :py:obj:`spatialetl.coverage.io.CoverageReader.CoverageReader`


   .. py:attribute:: nodefile


   .. py:attribute:: data


   .. py:attribute:: x_size


   .. py:attribute:: y_size


   .. py:method:: is_regular_grid()


   .. py:method:: get_x_size()


   .. py:method:: get_y_size()


   .. py:method:: get_t_size()


   .. py:method:: read_axis_x(xmin, xmax, ymin, ymax)


   .. py:method:: read_axis_y(xmin, xmax, ymin, ymax)


   .. py:method:: read_variable_bathymetry(xmin, xmax, ymin, ymax)


   .. py:method:: iocheck(fname)
      :classmethod:



   .. py:method:: read_swantable(fname, headers=[])

      Use this function to read data generated with the command table.
      Both NOHEAD and HEAD options can be read here.
      If using NOHEAD,
      the user must specify with variables are being read, for example:
      R  = swantools.io.SwanIO()
      headers = ["TIME","HSIGN","HSWELL","PDIR","DIR","TPS","PER","WINDX","WINDY","PROPAGAT"]
      df = R.read_swantable('file.txt',headers=headers)
      If usind HEAD option, just do:
      R  = swantools.io.SwanIO()
      df = R.read_swantable('file_with_headers.txt')
      The function will return a pandas DataFrame.



   .. py:method:: read_swanspc(fname)

      Function to read data generated with the SPECOUT command.
      The sixtase MUST be :
      'SPECOUT 'Location' SPEC2D ABS 'name.spc'
      Read the documentation in http://swanmodel.sourceforge.net to more details on spectral output.
      Inputs:
          fname:    the name of the file
      returns a xray DataArray



   .. py:method:: read_swanblock(fname, variables, stat=False)

      Function to read SWAN's BLOCK output statement.
      Both stationary and non-stationary versions can be handle.
      inputs :
          fname     = full name of the block file
          variables = list of variables to be read
          stat      = Whether stationary block or nor, default is False
      Obs: Curvilinear grids can't be read yet.
      Returns a xray DataSet



   .. py:method:: write_spectrum(fname, lat, lon, times, freqs, dirs, facs, spc)


   .. py:method:: write_tpar(data, fname='tpar')

      Given an array Nx5, where N is the number of timestamps and
      the columns are respecvely, TIME, HS, TP, DP, SP, this function
      will write TPAR formated  file named *fname*.



.. py:class:: Converters

   .. py:method:: np2nc(fname, lat, lon, ts, z, var)


   .. py:method:: spc2nc(fname, lat, lon, freq, dirr, time, facs, spc)


.. py:function:: swantime2datetime(time, inverse=False)

   Translating Swans's time strings to datetimes and vice-versa.
   See datetime module more information.


.. py:function:: dir2cat(theta)

   Given a array of directions, will return the respective
   categoreis (N,S,E,W, ect..). Credits to Eric Nardi.


.. py:function:: deg2uv(direction, intensity=False)

   Givean an array of directions will return the U and V
   components. intensity is optional.


