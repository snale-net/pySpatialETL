spatialetl.utils.bing.core.AerialImageRetrieval
===============================================

.. py:module:: spatialetl.utils.bing.core.AerialImageRetrieval

.. autoapi-nested-parse::

   __author__ = Linlin Chen
   __email__ = lchen96@hawk.iit.edu


   @Description:
   This module is used to retrieve satellite/aerial image.
   Given a bounding box, which is composed of left up corner coordinate (latitude, longitude)
   and right down corner coordinate (latitude, longitude).
   Return an aerial imagery (with maximum resolution available) downloaded from Bing map tile system.

   @Edit:
   20/11/2018 Refactoring by Fabien Rétif - fabien.retif@zoho.com



Classes
-------

.. autoapisummary::

   spatialetl.utils.bing.core.AerialImageRetrieval.AerialImageRetrieval


Module Contents
---------------

.. py:class:: AerialImageRetrieval(lat1, lon1, lat2, lon2, outputfile)

   Bases: :py:obj:`object`


   The class for aerial image retrieval

   To create an AerialImageRetrieval object, simply give upper left latitude, longitude,
   and lower right latitude and longitude


   .. py:attribute:: BASEURL
      :value: 'http://h0.ortho.tiles.virtualearth.net/tiles/a{0}.jpeg?g=131'



   .. py:attribute:: IMAGEMAXSIZE
      :value: 11366912



   .. py:attribute:: TILESIZE
      :value: 256



   .. py:attribute:: lat1


   .. py:attribute:: lon1


   .. py:attribute:: lat2


   .. py:attribute:: lon2


   .. py:attribute:: outputFilename


   .. py:attribute:: tempFile


   .. py:method:: download_image(quadkey)

      This method is used to download a tile image given the quadkey from Bing tile system

      :param quadkey {[string]} -- [The quadkey for a tile image]:

      :returns: [Image] -- [A PIL Image]



   .. py:method:: is_valid_image(image)

      Check whether the downloaded image is valid,
      by comparing the downloaded image with a NULL image returned by any unsuccessfully retrieval

      Bing tile system will return the same NULL image if the query quadkey is not existed in the Bing map database.

      :param image {[Image]} -- [a Image type image to be valided]:

      :returns: [boolean] -- [whether the image is valid]



   .. py:method:: max_resolution_imagery_retrieval()

      The main aerial retrieval method

      It will firstly determine the appropriate level used to retrieve the image.
      The appropriate level should satisfy:
          1. All the tile image within the given bounding box at that level should all exist
          2. The retrieved image cannot exceed the maximum supported image size, which is 8192*8192 (Otherwise the image size will be too large if the bounding box is very large)

      Then for the given level, we can download each aerial tile image, and stitch them together.

      Lastly, we have to crop the image based on the given bounding box

      :returns: [boolean] -- [indicate whether the aerial image retrieval is successful]



   .. py:method:: horizontal_retrieval_and_stitch_image(tileX_start, tileX_end, tileY, level)

      Horizontally retrieve tile images and then stitch them together,
      start from tileX_start and end at tileX_end, tileY will remain the same

      :param tileX_start {[int]} -- [the starting tileX index]:
      :param tileX_end {[int]} -- [the ending tileX index]:
      :param tileY {[int]} -- [the tileY index]:
      :param level {[int]} -- [level used to retrieve image]:

      :returns: [boolean, Image] -- [whether such retrieval is successful; If successful, returning the stitched image, otherwise None]



