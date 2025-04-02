spatialetl.utils.bing.core.ogr2ogr
==================================

.. py:module:: spatialetl.utils.bing.core.ogr2ogr


Attributes
----------

.. autoapisummary::

   spatialetl.utils.bing.core.ogr2ogr.nLastTick
   spatialetl.utils.bing.core.ogr2ogr.bSkipFailures
   spatialetl.utils.bing.core.ogr2ogr.nGroupTransactions
   spatialetl.utils.bing.core.ogr2ogr.bPreserveFID
   spatialetl.utils.bing.core.ogr2ogr.nFIDToFetch
   spatialetl.utils.bing.core.ogr2ogr.GeomOperation
   spatialetl.utils.bing.core.ogr2ogr.version_num


Classes
-------

.. autoapisummary::

   spatialetl.utils.bing.core.ogr2ogr.ScaledProgressObject
   spatialetl.utils.bing.core.ogr2ogr.TargetLayerInfo
   spatialetl.utils.bing.core.ogr2ogr.AssociatedLayers
   spatialetl.utils.bing.core.ogr2ogr.Enum


Functions
---------

.. autoapisummary::

   spatialetl.utils.bing.core.ogr2ogr.ScaledProgressFunc
   spatialetl.utils.bing.core.ogr2ogr.EQUAL
   spatialetl.utils.bing.core.ogr2ogr.TermProgress
   spatialetl.utils.bing.core.ogr2ogr.main
   spatialetl.utils.bing.core.ogr2ogr.Usage
   spatialetl.utils.bing.core.ogr2ogr.CSLFindString
   spatialetl.utils.bing.core.ogr2ogr.IsNumber
   spatialetl.utils.bing.core.ogr2ogr.LoadGeometry
   spatialetl.utils.bing.core.ogr2ogr.wkbFlatten
   spatialetl.utils.bing.core.ogr2ogr.SetZ
   spatialetl.utils.bing.core.ogr2ogr.SetupTargetLayer
   spatialetl.utils.bing.core.ogr2ogr.TranslateLayer


Module Contents
---------------

.. py:class:: ScaledProgressObject(min, max, cbk, cbk_data=None)

   .. py:attribute:: min


   .. py:attribute:: max


   .. py:attribute:: cbk


   .. py:attribute:: cbk_data
      :value: None



.. py:function:: ScaledProgressFunc(pct, msg, data)

.. py:function:: EQUAL(a, b)

.. py:data:: nLastTick
   :value: -1


.. py:function:: TermProgress(dfComplete, pszMessage, pProgressArg)

.. py:class:: TargetLayerInfo

   .. py:attribute:: poDstLayer
      :value: None



   .. py:attribute:: poCT
      :value: None



   .. py:attribute:: panMap
      :value: None



   .. py:attribute:: iSrcZField
      :value: None



.. py:class:: AssociatedLayers

   .. py:attribute:: poSrcLayer
      :value: None



   .. py:attribute:: psInfo
      :value: None



.. py:data:: bSkipFailures
   :value: False


.. py:data:: nGroupTransactions
   :value: 200


.. py:data:: bPreserveFID
   :value: False


.. py:data:: nFIDToFetch

.. py:class:: Enum

   Bases: :py:obj:`set`


   set() -> new empty set object
   set(iterable) -> new set object

   Build an unordered collection of unique elements.


   .. py:method:: __getattr__(name)


.. py:data:: GeomOperation

.. py:function:: main(args=None, progress_func=TermProgress, progress_data=None)

.. py:function:: Usage()

.. py:function:: CSLFindString(v, mystr)

.. py:function:: IsNumber(pszStr)

.. py:function:: LoadGeometry(pszDS, pszSQL, pszLyr, pszWhere)

.. py:function:: wkbFlatten(x)

.. py:function:: SetZ(poGeom, dfZ)

.. py:function:: SetupTargetLayer(poSrcDS, poSrcLayer, poDstDS, papszLCO, pszNewLayerName, bTransform, poOutputSRS, bNullifyOutputSRS, poSourceSRS, papszSelFields, bAppend, eGType, bPromoteToMulti, nCoordDim, bOverwrite, papszFieldTypesToString, bWrapDateline, bExplodeCollections, pszZField, pszWHERE)

.. py:function:: TranslateLayer(psInfo, poSrcDS, poSrcLayer, poDstDS, poOutputSRS, bNullifyOutputSRS, eGType, bPromoteToMulti, nCoordDim, eGeomOp, dfGeomOpParam, nCountLayerFeatures, poClipSrc, poClipDst, bExplodeCollections, nSrcFileSize, pnReadFeatureCount, pfnProgress, pProgressArg)

.. py:data:: version_num

