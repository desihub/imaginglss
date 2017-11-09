.. ImagingLSS documentation master file, created by
   sphinx-quickstart on Sun Mar 15 12:45:33 2015.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

Welcome to ImagingLSS's documentation!
======================================

ImagingLSS implements a data model of the DECALS imaging survey data.
The intended use of ImagingLSS is to generate large scale structure 
catalogues and masks from the DECALS imaging survey data.

Conceptually, ImagingLSS is a consumer of the products of `Legacy Pipe <https://github.com/legacysurvey/legacypipe>`_, 
the data pipe line of the DECam Legacy Survey.

Currently, `DR2 <https://desi.lbl.gov/trac/wiki/DecamLegacy/DR2>`_ is supported.
Previous Data Releases are also supported (since Early Data Release, or EDR), 
though in the future we may drop the support in order to simplify the code base.

Users Guide
===========
.. toctree::
   :maxdepth: 2

   install
   datapipeline
   datamodel
   helpers
   examples

API Reference
=============
.. toctree::
   :maxdepth: 2

   modules

Known Issues
============

- Due to internal caching of the catalogue, the running time is faster after first time.

- If the data release has updated (rare but it happened), 
  and the cache files are out of date. 
  Run forget_cache.py to clear the cache.

Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`


