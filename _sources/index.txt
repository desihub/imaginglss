.. ImagingLSS documentation master file, created by
   sphinx-quickstart on Sun Mar 15 12:45:33 2015.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

Welcome to ImagingLSS's documentation!
======================================

This is a data model of the DECAM imaging survey data 
and potentially also for other desi targeting imaging survey, 
for Large scale structure analysis.

The modelled data is a product from the Tractor, the source identification
software used in by DESI.

We have support on 
`Early Data Release (EDR) data <https://desi.lbl.gov/trac/wiki/DecamLegacy/EDRfiles>`_
and `DR1 <https://desi.lbl.gov/trac/wiki/DecamLegacy/DR1>`_ .

Users Guide
===========
.. toctree::
   :maxdepth: 2

   install
   datamodel
   helpers
   examples
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


