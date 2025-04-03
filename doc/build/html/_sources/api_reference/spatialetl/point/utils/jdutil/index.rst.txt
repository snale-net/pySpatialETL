spatialetl.point.utils.jdutil
=============================

.. py:module:: spatialetl.point.utils.jdutil

.. autoapi-nested-parse::

   Functions for converting dates to/from JD and MJD. Assumes dates are historical
   dates, including the transition from the Julian calendar to the Gregorian
   calendar in 1582. No support for proleptic Gregorian/Julian calendars.

   :Author: Matt Davis
   :Website: http://github.com/jiffyclub



Classes
-------

.. autoapisummary::

   spatialetl.point.utils.jdutil.datetime


Functions
---------

.. autoapisummary::

   spatialetl.point.utils.jdutil.mjd_to_jd
   spatialetl.point.utils.jdutil.jd_to_mjd
   spatialetl.point.utils.jdutil.date_to_jd
   spatialetl.point.utils.jdutil.jd_to_date
   spatialetl.point.utils.jdutil.hmsm_to_days
   spatialetl.point.utils.jdutil.days_to_hmsm
   spatialetl.point.utils.jdutil.datetime_to_jd
   spatialetl.point.utils.jdutil.jd_to_datetime
   spatialetl.point.utils.jdutil.timedelta_to_days


Module Contents
---------------

.. py:function:: mjd_to_jd(mjd)

   Convert Modified Julian Day to Julian Day.

   :param mjd: Modified Julian Day
   :type mjd: float

   :returns: **jd** -- Julian Day
   :rtype: float


.. py:function:: jd_to_mjd(jd)

   Convert Julian Day to Modified Julian Day

   :param jd: Julian Day
   :type jd: float

   :returns: **mjd** -- Modified Julian Day
   :rtype: float


.. py:function:: date_to_jd(year, month, day)

   Convert a date to Julian Day.

   Algorithm from 'Practical Astronomy with your Calculator or Spreadsheet',
       4th ed., Duffet-Smith and Zwart, 2011.

   :param year: Year as integer. Years preceding 1 A.D. should be 0 or negative.
                The year before 1 A.D. is 0, 10 B.C. is year -9.
   :type year: int
   :param month: Month as integer, Jan = 1, Feb. = 2, etc.
   :type month: int
   :param day: Day, may contain fractional part.
   :type day: float

   :returns: **jd** -- Julian Day
   :rtype: float

   .. rubric:: Examples

   Convert 6 a.m., February 17, 1985 to Julian Day

   >>> date_to_jd(1985,2,17.25)
   2446113.75


.. py:function:: jd_to_date(jd)

   Convert Julian Day to date.

   Algorithm from 'Practical Astronomy with your Calculator or Spreadsheet',
       4th ed., Duffet-Smith and Zwart, 2011.

   :param jd: Julian Day
   :type jd: float

   :returns: * **year** (*int*) -- Year as integer. Years preceding 1 A.D. should be 0 or negative.
               The year before 1 A.D. is 0, 10 B.C. is year -9.
             * **month** (*int*) -- Month as integer, Jan = 1, Feb. = 2, etc.
             * **day** (*float*) -- Day, may contain fractional part.

   .. rubric:: Examples

   Convert Julian Day 2446113.75 to year, month, and day.

   >>> jd_to_date(2446113.75)
   (1985, 2, 17.25)


.. py:function:: hmsm_to_days(hour=0, min=0, sec=0, micro=0)

   Convert hours, minutes, seconds, and microseconds to fractional days.

   :param hour: Hour number. Defaults to 0.
   :type hour: int, optional
   :param min: Minute number. Defaults to 0.
   :type min: int, optional
   :param sec: Second number. Defaults to 0.
   :type sec: int, optional
   :param micro: Microsecond number. Defaults to 0.
   :type micro: int, optional

   :returns: **days** -- Fractional days.
   :rtype: float

   .. rubric:: Examples

   >>> hmsm_to_days(hour=6)
   0.25


.. py:function:: days_to_hmsm(days)

   Convert fractional days to hours, minutes, seconds, and microseconds.
   Precision beyond microseconds is rounded to the nearest microsecond.

   :param days: A fractional number of days. Must be less than 1.
   :type days: float

   :returns: * **hour** (*int*) -- Hour number.
             * **min** (*int*) -- Minute number.
             * **sec** (*int*) -- Second number.
             * **micro** (*int*) -- Microsecond number.

   :raises ValueError: If `days` is >= 1.

   .. rubric:: Examples

   >>> days_to_hmsm(0.1)
   (2, 24, 0, 0)


.. py:function:: datetime_to_jd(date)

   Convert a `datetime.datetime` object to Julian Day.

   :param date:
   :type date: `datetime.datetime` instance

   :returns: **jd** -- Julian day.
   :rtype: float

   .. rubric:: Examples

   >>> d = datetime.datetime(1985,2,17,6)
   >>> d
   datetime.datetime(1985, 2, 17, 6, 0)
   >>> jdutil.datetime_to_jd(d)
   2446113.75


.. py:function:: jd_to_datetime(jd)

   Convert a Julian Day to an `jdutil.datetime` object.

   :param jd: Julian day.
   :type jd: float

   :returns: **dt** -- `jdutil.datetime` equivalent of Julian day.
   :rtype: `jdutil.datetime` object

   .. rubric:: Examples

   >>> jd_to_datetime(2446113.75)
   datetime(1985, 2, 17, 6, 0)


.. py:function:: timedelta_to_days(td)

   Convert a `datetime.timedelta` object to a total number of days.

   :param td:
   :type td: `datetime.timedelta` instance

   :returns: **days** -- Total number of days in the `datetime.timedelta` object.
   :rtype: float

   .. rubric:: Examples

   >>> td = datetime.timedelta(4.5)
   >>> td
   datetime.timedelta(4, 43200)
   >>> timedelta_to_days(td)
   4.5


.. py:class:: datetime

   Bases: :py:obj:`datetime.datetime`


   A subclass of `datetime.datetime` that performs math operations by first
   converting to Julian Day, then back to a `jdutil.datetime` object.

   Addition works with `datetime.timedelta` objects, subtraction works with
   `datetime.timedelta`, `datetime.datetime`, and `jdutil.datetime` objects.
   Not all combinations work in all directions, e.g.
   `timedelta - datetime` is meaningless.

   .. seealso::

      :py:obj:`datetime.datetime`
          Parent class.


   .. py:method:: __add__(other)

      Add a datetime and a timedelta.



   .. py:method:: __radd__(other)


   .. py:method:: __sub__(other)

      Subtract two datetimes, or a datetime and a timedelta.



   .. py:method:: __rsub__(other)


   .. py:method:: to_jd()

      Return the date converted to Julian Day.




   .. py:method:: to_mjd()

      Return the date converted to Modified Julian Day.




