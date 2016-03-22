#!/usr/bin/env python
#
# Code to query the VETO mask of objects/randoms
# It takes the NOISES extension as an input
# It writers a VETO extension.

# Usage, see python query_veto.py -h
#
from __future__ import print_function


__author__ = "Yu Feng and Martin White"
__version__ = "1.0"
__email__  = "yfeng1@berkeley.edu or mjwhite@lbl.gov"

import numpy as np
from imaginglss             import DECALS
from imaginglss.analysis    import tycho_veto
from imaginglss.utils       import output

from argparse import ArgumentParser

ap = ArgumentParser("select_objs.py")
ap.add_argument("input", type=output.writer)
ap.add_argument("output", type=output.writer)
ap.add_argument("--conf", default=None,
        help="Path to the imaginglss config file, default is from DECALS_PY_CONFIG")

ns = ap.parse_args()
ns.conf = DECALS(ns.conf)

from mpi4py import MPI

np.seterr(divide='ignore', invalid='ignore')

def query_veto(ns, comm=MPI.COMM_WORLD):
    """

    """
    #FIXME: this reads only text files.
    objects = ns.input.read('NOISES')
    RA = objects['RA']
    DEC = objects['DEC']
    allvetos = [i for i in dir(tycho_veto) if not str(i).startswith( '_' )]
    dataset = np.zeros(len(RA), dtype=
            [('VETO', 'uint64')] +
            list(zip(allvetos, ['?'] * len(allvetos)))
            )

    for vetoname in allvetos:
        veto = getattr(tycho_veto, vetoname)
        mask = veto(ns.conf.tycho, (RA, DEC))
        dataset['VETO'][mask] |= veto.BIT
        dataset[vetoname][mask] = True

    return dataset

if __name__=="__main__":

    VETO = query_veto(ns)

    if MPI.COMM_WORLD.rank == 0:
        ns.output.write(VETO, ns.__dict__, 'VETO')
