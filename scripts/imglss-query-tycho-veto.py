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
import h5py

from imaginglss             import DECALS
from imaginglss.analysis    import tycho_veto
from imaginglss.analysis    import veto

from imaginglss.cli import CLI

cli = CLI(
"""
Query the TYCHOVETO flags of input data. The position is taken from the NOISES extension of input.
The result is written to the TYCHOVETO extension of output.

Currently, only veto by proximity to tycho stars are implemented. Each veto in
imaginglss.analysis.tycho_veto is calculated as a column in the TYCHOVETO extension.

Unfortunately, this script is not sufficiently smart to decide the correct TYCHOVETO for the target type.
Therefore, no combined veto flag is generated.
"""
)
cli.add_argument("catalogue", help="HDF5 catalogue file, can be either random or objects. TYCHO_VETO dataset will be added ")

ns = cli.parse_args()
decals = DECALS(ns.conf)

np.seterr(divide='ignore', invalid='ignore')

def query_veto(decals, ns):
    """
        calculate VETO flag for all proximity vetos defined in tycho_veto.
    """

    with h5py.File(ns.catalogue, 'r') as ff:
        RA = ff['RA'][:]
        DEC = ff['DEC'][:]

    allvetos = [i for i in dir(tycho_veto) if not str(i).startswith( '_' )]
    dataset = np.zeros(len(RA), dtype=
            list(zip(allvetos, ['?'] * len(allvetos)))
            )

    for ibit, vetoname in enumerate(allvetos):
        vetotype = getattr(tycho_veto, vetoname)
        R = vetotype(decals.tycho)
        centers = (decals.tycho['RA'], decals.tycho['DEC'])
        mask = veto.veto((RA, DEC), centers, R)
        # if we want to combine the bits, do it here.
        # but there is no point of doing so for all tycho based proximity
        # vetos. we will assembly the full veto bitmask later in the pipeline.
        dataset[vetoname][mask] = True
        print(vetoname, dataset[vetoname].sum())
    return dataset

if __name__=="__main__":

    VETO = query_veto(decals, ns)

    with h5py.File(ns.catalogue, 'r+') as ff:
        if 'TYCHO_VETO' in ff:
            del ff['TYCHO_VETO']
        ds = ff.create_dataset('TYCHO_VETO', data=VETO)
