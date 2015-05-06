#!/usr/bin/env python
#
# Code to build the catalogue cache
#
# Usage: python build_cache.py
#
from __future__ import print_function
from sys import stdout

__author__ = "Yu Feng and Martin White"
__version__ = "1.0"
__email__  = "yfeng1@berkeley.edu or mjwhite@lbl.gov"

import os.path; import sys; sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..'))

from imaginglss             import DECALS

decals = DECALS()
def report(processed, total):
    stdout.write("processing %d/%d ... \n" %( processed, total))
    stdout.flush()
decals.datarelease.catalogue.build_cache(report=report)
print("done")

