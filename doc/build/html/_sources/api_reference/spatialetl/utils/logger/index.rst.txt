spatialetl.utils.logger
=======================

.. py:module:: spatialetl.utils.logger


Attributes
----------

.. autoapisummary::

   spatialetl.utils.logger.logging


Classes
-------

.. autoapisummary::

   spatialetl.utils.logger.RunFilter
   spatialetl.utils.logger.myLogger


Module Contents
---------------

.. py:class:: RunFilter(name='')

   Bases: :py:obj:`logging.Filter`


   Filter instances are used to perform arbitrary filtering of LogRecords.

   Loggers and Handlers can optionally use Filter instances to filter
   records as desired. The base filter class only allows events which are
   below a certain point in the logger hierarchy. For example, a filter
   initialized with "A.B" will allow events logged by loggers "A.B",
   "A.B.C", "A.B.C.D", "A.B.D" etc. but not "A.BB", "B.A.B" etc. If
   initialized with the empty string, all events are passed.


   .. py:method:: filter(record)

      Determine if the specified record is to be logged.

      Returns True if the record should be logged, or False otherwise.
      If deemed appropriate, the record may be modified in-place.



.. py:class:: myLogger(name)

   Bases: :py:obj:`logging.Logger`


   Instances of the Logger class represent a single logging channel. A
   "logging channel" indicates an area of an application. Exactly how an
   "area" is defined is up to the application developer. Since an
   application can have any number of areas, logging channels are identified
   by a unique string. Application areas can be nested (e.g. an area
   of "input processing" might include sub-areas "read CSV files", "read
   XLS files" and "read Gnumeric files"). To cater for this natural nesting,
   channel names are organized into a namespace hierarchy where levels are
   separated by periods, much like the Java or Python package namespace. So
   in the instance given above, channel names might be "input" for the upper
   level, and "input.csv", "input.xls" and "input.gnu" for the sub-levels.
   There is no arbitrary limit to the depth of nesting.


   .. py:attribute:: TIMING
      :value: 200



   .. py:attribute:: RUN
      :value: 25



   .. py:attribute:: INFO
      :value: 20



   .. py:attribute:: DEBUG
      :value: 10



   .. py:attribute:: WARNING
      :value: 30



   .. py:attribute:: ERROR
      :value: 40



   .. py:attribute:: CRITICAL
      :value: 50



   .. py:method:: debug(msg, *args, **kwargs)

      Log 'msg % args' with severity 'DEBUG'.

      To pass exception information, use the keyword argument exc_info with
      a true value, e.g.

      logger.debug("Houston, we have a %s", "thorny problem", exc_info=True)



   .. py:method:: info(msg, *args, **kwargs)

      Log 'msg % args' with severity 'INFO'.

      To pass exception information, use the keyword argument exc_info with
      a true value, e.g.

      logger.info("Houston, we have a %s", "notable problem", exc_info=True)



   .. py:method:: warning(msg, *args, **kwargs)

      Log 'msg % args' with severity 'WARNING'.

      To pass exception information, use the keyword argument exc_info with
      a true value, e.g.

      logger.warning("Houston, we have a %s", "bit of a problem", exc_info=True)



   .. py:method:: timing(msg, *args, **kwargs)


   .. py:method:: setLevel(level)

      Set the logging level of this logger.  level must be an int or a str.



.. py:data:: logging

