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
from imaginglss.utils import filehandler
from mpi4py import MPI
import numpy

comm = MPI.COMM_WORLD

decals = DECALS()
cat = decals.datarelease.catalogue

filenames = cat.filenames

mystart = len(filenames) * comm.rank // comm.size
myend = len(filenames) * (comm.rank + 1)// comm.size

if comm.rank == 0:
    print("reading files ...")
data = cat.build_cache(filenames[mystart:myend])

comm.barrier()

if comm.rank == 0:
    print("writing cache ...")
for i in range(comm.size):
    comm.barrier()
    if i != comm.rank: continue
    mode = 'w' if comm.rank == 0 else 'a'
    filehandler.write(cat.cachedir, data, mode=mode)
    print("chunk %d / %d done" %( comm.rank, comm.size))

comm.barrier()

if comm.rank == 0:
    total = len(filenames)
    d = {
        'nfiles': numpy.array([total],dtype='i8') 
        }

    filehandler.write(cat.cachedir, d)

