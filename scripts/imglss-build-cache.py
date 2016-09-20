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

from imaginglss import DECALS
import numpy

from imaginglss.cli import CLI
from imaginglss.analysis import cache
from imaginglss.model.catalogue import Catalogue

ap = CLI("Build cache")

ns = ap.parse_args()

decals = DECALS(ns.conf)
builder = cache.CacheBuilder(decals.sweep_dir, decals.cache_dir, Catalogue.COLUMNS)

builder.build()

