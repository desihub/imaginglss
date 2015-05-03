#!/usr/bin/env python
#
# Python script to generate clear the internally cached catalogue data
# these files are stored in $DECALS_CACHE directory.
#

from __future__ import print_function

__author__ = "Yu Feng and Martin White"
__version__ = "1.0"
__email__  = "yfeng1@berkeley.edu or mwhite@berkeley.edu"

import os.path; import sys; sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..'))

from imaginglss.datarelease import DataRelease

dr = DataRelease()
cat = dr.catalogue

for key in cat:
    cat.forget(key)

