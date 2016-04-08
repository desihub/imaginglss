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
from imaginglss.utils import filehandler
import numpy

from argparse import ArgumentParser

ap = ArgumentParser()
ap.add_argument("--conf", default=None, 
        help="Path to the imaginglss config file, default is from DECALS_PY_CONFIG")

ns = ap.parse_args()

from mpi4py import MPI
comm = MPI.COMM_WORLD

# Staggered start up builds covered_bricks.i8 on root rank only
if comm.rank != 0:
    comm.barrier()

decals = DECALS(ns.conf)
cat = decals.datarelease.catalogue

# Staggered start up builds covered_bricks.i8 on root rank only
if comm.rank == 0:
    comm.barrier()

filenames = list(cat.filenames.values())

mystart = len(filenames) * comm.rank // comm.size
myend = len(filenames) * (comm.rank + 1)// comm.size

if comm.rank == 0:
    print("reading files ...")

data = cat.build_cache(filenames[mystart:myend])

comm.barrier()

if comm.rank == 0:
    print("writing cache ...")

if comm.rank == 0:
    filehandler.write(cat.cachedir, data, mode='w')

comm.barrier()

N = comm.allgather(len(data))

filehandler.write(cat.cachedir, data, mode='r+', offset=sum(N[:comm.rank]))


print("chunk %d / %d done" %( comm.rank, comm.size))
comm.barrier()

if comm.rank == 0:
    print("%d items are written" % sum(N))
    total = len(filenames)
    d = {
        'nfiles': numpy.array([total],dtype='i8') 
        }

    filehandler.write(cat.cachedir, d)

