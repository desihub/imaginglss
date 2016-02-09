Helper routines
===============

.. contents::
    :depth: 2
 
In order to reduce the number of dependecies, 
we implemented a few coordinate transformations of the WCS system. Only
relevant special cases are implemented.

There is also a sub-package for in-node parallel multiprocess. 

WCS
---

 - :py:mod:`imaginglss.utils.wcs_tangent` impments the WCS tangent format.
 - :py:mod:`imaginglss.utils.wcs_simplezea` implement the WCS ZEA format; oriented at north / south pole
 - :py:mod:`imaginglss.utils.euler` implements the transformation between Galactic and RA/DEC. coordinates.

Sharedmem
---------

 - :py:mod:`imaginglss.utils.sharedmem` implements a MapReduce pool object for in-node parallellism. See :doc:`sharedmem`

