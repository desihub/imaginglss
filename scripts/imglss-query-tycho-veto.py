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

ap = ArgumentParser(
description=
"""
Query the TYCHOVETO flags of input data. The position is taken from the NOISES extension of input.
The result is written to the TYCHOVETO extension of output.

Currently, only veto by proximity to tycho stars are implemented. Each veto in
imaginglss.analysis.tycho_veto is calculated as a column in the TYCHOVETO extension.

Unfortunately, this script is not sufficiently smart to decide the correct TYCHOVETO for the target type.
Therefore, no combined veto flag is generated.
"""
)
ap.add_argument("input", type=output.writer, 
    help="Reads the position from the NOISES extension.")
ap.add_argument("output", type=output.writer,
    help="Writes the veto flags to the TYCHOVETO extension")
ap.add_argument("--conf", default=None,
        help="Path to the imaginglss config file, default is from DECALS_PY_CONFIG")

ns = ap.parse_args()
ns.conf = DECALS(ns.conf)

from mpi4py import MPI

np.seterr(divide='ignore', invalid='ignore')

def query_veto(ns, comm=MPI.COMM_WORLD):
    """
        calculate VETO flag for all proximity vetos defined in tycho_veto.
    """

    objects = ns.input.read('NOISES')
    RA = objects['RA']
    DEC = objects['DEC']
    allvetos = [i for i in dir(tycho_veto) if not str(i).startswith( '_' )]
    dataset = np.zeros(len(RA), dtype=
            list(zip(allvetos, ['?'] * len(allvetos)))
            )

    for ibit, vetoname in enumerate(allvetos):
        veto = getattr(tycho_veto, vetoname)
        mask = ns.conf.tycho.veto((RA, DEC), veto)
        # if we want to combine the bits, do it here.
        # but there is no point of doing so for all tycho based proximity
        # vetos. we will assembly the full veto bitmask later in the pipeline.
        dataset[vetoname][mask] = True

    return dataset

if __name__=="__main__":

    VETO = query_veto(ns)

    if MPI.COMM_WORLD.rank == 0:
        ns.output.write(VETO, ns.__dict__, 'TYCHOVETO')
