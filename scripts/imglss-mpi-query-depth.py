#!/usr/bin/env python

from __future__ import print_function,division

__author__ = "Yu Feng and Martin White"
__version__ = "1.0"
__email__  = "yfeng1@berkeley.edu or mjwhite@lbl.gov"

import h5py

from   imaginglss             import DECALS

from imaginglss.cli import CLI

cli = CLI("""
Query Depth from DECALS data for input RA DEC of points.
The input must be saved in a HDF5 with two datasets 'RA' and 'DEC'.
The output will be written in the same file as INTRINSIC_NOISELEVEL data set.
To lookup the columns, use the dictionary in `imaginglss.model.datcliroduct.bands`.

The output of this script can be directly fed into imglss-query-completeness.py
as the query input.

""")

cli.add_argument("query", help="An HDF5 file with RA and DEC dataset, the position of to query the depth." )

ns = cli.parse_args()
decals = DECALS(ns.conf)

import numpy             as np
from   imaginglss.model.datarelease import Footprint
from   imaginglss.model import dataproduct
from   mpi4py            import MPI
from   imaginglss.analysis    import cuts

np.seterr(divide='ignore', invalid='ignore')


def query_depth(decals, ns, comm=MPI.COMM_WORLD):
    """
        query the depth.
    """
    with h5py.File(ns.query, 'r') as ff:
        RA = ff['RA']
        DEC = ff['DEC']
        mystart = RA.size * comm.rank // comm.size
        myend = RA.size * (comm.rank + 1) // comm.size
        coord = (ff['RA'][mystart:myend], ff['DEC'][mystart:myend])

    # Get the total footprint bounds, to throw randoms within, and an E(B-V)
    # map instance.
    dr = decals.datarelease
    sfd= decals.sfdmap

    randoms = np.empty(len(coord[0]), dtype=dataproduct.RandomCatalogue)
    randoms['RA'] = coord[0]
    randoms['DEC'] = coord[1]

    ns.Nran = sum(comm.allgather(len(coord[0])))

    cat_lim = dr.read_depths(coord, 'grz')

    # It's also useful to the 1 sigma limits later.
    randoms['INTRINSIC_NOISELEVEL'][:, :6] = (cat_lim['DECAM_DEPTH'] ** -0.5 / cat_lim['DECAM_MW_TRANSMISSION'])
    randoms['INTRINSIC_NOISELEVEL'][:, 6:] = 0

    nanmask = np.isnan(randoms['INTRINSIC_NOISELEVEL'])
    randoms['INTRINSIC_NOISELEVEL'][nanmask] = np.inf

    randoms = comm.gather(randoms)

    if comm.rank == 0:
        randoms = np.concatenate(randoms)
        print(len(randoms), 'Done!')
    return randoms

if __name__ == '__main__':

    randoms = query_depth(decals, ns)

    if MPI.COMM_WORLD.rank == 0:
        with h5py.File(ns.query, 'r+') as ff:
            if 'INTRINSIC_NOISELEVEL' in ff:
                del ff['INTRINSIC_NOISELEVEL']
            ds = ff.create_dataset('INTRINSIC_NOISELEVEL', data=randoms['INTRINSIC_NOISELEVEL'])
            ds.attrs.update(cli.prune_namespace(ns))
