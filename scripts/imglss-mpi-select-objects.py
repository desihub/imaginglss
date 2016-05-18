#!/usr/bin/env python
#
# Code to select targets from Tractor catalogs.
# The target selection criteria are described at:
#   https://desi.lbl.gov/trac/wiki/TargetSelection
#
# These catalogs also need to be vetod using the
# depths and masks which are part of Tractor, since
# catalog entries can appear in places where the
# nominal depth is insufficient to be complete.
#
# Usage: python select_objs.py [--plot]
#
from __future__ import print_function


__author__ = "Yu Feng and Martin White"
__version__ = "1.0"
__email__  = "yfeng1@berkeley.edu or mjwhite@lbl.gov"

import numpy as np
import h5py

from imaginglss             import DECALS
from imaginglss.analysis    import targetselection
from imaginglss.analysis    import cuts
from imaginglss.model       import dataproduct
from imaginglss.cli         import CLI

cli = CLI("Select Objects based on Target definitions", enable_target_plugins=True)
cli.add_argument("--use-tractor-depth", action='store_true', default=False, help="Use Tractor's Depth in the catalogue, very fast!")
cli.add_argument("ObjectType", choices=[i for i in targetselection.__all__])
cli.add_argument("output", help="Output file name. A new object catalogue file will be created.")

ns = cli.parse_args()

decals = DECALS(ns.conf)

from mpi4py import MPI

np.seterr(divide='ignore', invalid='ignore')

def select_objs(decals, ns, comm=MPI.COMM_WORLD):
    """
    Does the actual selection, imposing cuts on the fluxes
    """

    dtype1 = np.dtype([
        ('RA', 'f8'),
        ('DEC', 'f8'),
        ('PHOTO_Z', 'f8'),
        ('DECAM_INTRINSIC_FLUX', ('f4', 6))])

    dtype2 = np.dtype([
        ('RA', 'f8'),
        ('DEC', 'f8'),
        ('DECAM_INTRINSIC_NOISE_LEVEL', ('f4', 6)),
        ])

    dtype3 = np.dtype([
        ('DECAM_CONFIDENCE', ('f4', 6)),
        ])

    # Get instances of a data release and SFD dust map.
    dr     = decals.datarelease
    sfd    = decals.sfdmap
    cat    = dr.catalogue
    #
    mystart = cat.size * comm.rank // comm.size
    myend = cat.size * (comm.rank + 1) // comm.size
    #
    mine = slice(mystart, myend)

    with dr.catalogue as cat:
        rows = cat[mine]
        fluxcut = getattr(targetselection, ns.ObjectType)
        mask = cuts.apply(comm, fluxcut, rows)

    if comm.rank == 0:
        print('Rank 0 with', myend - mystart, 'items')
        print('Rank 0 selected', mask.sum(), 'items')

    targets = np.empty(mask.sum(), dtype=dataproduct.ObjectCatalogue)

    targets['RA']   = cat[ 'RA'][mine][mask]
    targets['DEC']   = cat['DEC'][mine][mask]

    DECAM_FLUX = cat['DECAM_FLUX'][mine][mask] / cat['DECAM_MW_TRANSMISSION'][mine][mask]
    WISE_FLUX = cat['WISE_FLUX'][mine][mask] / cat['WISE_MW_TRANSMISSION'][mine][mask]
    targets['INTRINSIC_FLUX'][:, :6] = (DECAM_FLUX).clip(1e-15,1e15)
    targets['INTRINSIC_FLUX'][:, 6:] = (WISE_FLUX).clip(1e-15,1e15)

    # Now we need to pass this through our mask since galaxies can
    # appear even in regions where our nominal depth is insufficient
    # for a complete ns.ObjectTypee.

    # note that although DR2 has these numbers, we query the raw images
    # to keep it consistent with make_randoms.py.  This shall agree with
    # the DR2 tractor file numbers.

    # ... and extract only the objects which passed the cuts.

    cat_lim = np.empty(len(targets), dtype=[
        ('DECAM_DEPTH', cat.dtype['DECAM_DEPTH']),
        ('WISE_FLUX_IVAR', cat.dtype['WISE_FLUX_IVAR']),
        ('DECAM_MW_TRANSMISSION', cat.dtype['DECAM_MW_TRANSMISSION']),
        ('WISE_MW_TRANSMISSION', cat.dtype['WISE_MW_TRANSMISSION']),
        ])
    if ns.use_tractor_depth:
        cat_lim['DECAM_DEPTH'][:] = cat['DECAM_DEPTH'][mine][mask]
        cat_lim['WISE_FLUX_IVAR'][:] = cat['WISE_FLUX_IVAR'][mine][mask]
        cat_lim['DECAM_MW_TRANSMISSION'][:] = cat['DECAM_MW_TRANSMISSION'][mine][mask]
        cat_lim['WISE_MW_TRANSMISSION'][:] = cat['WISE_MW_TRANSMISSION'][mine][mask]
    else:
        cat_lim1 = dr.read_depths((targets['RA'], targets['DEC']), 'grz')
        cat_lim['DECAM_DEPTH'][:] = cat_lim1['DECAM_DEPTH']
        cat_lim['WISE_FLUX_IVAR'][:] = cat['WISE_FLUX_IVAR'][mine][mask]
        cat_lim['DECAM_MW_TRANSMISSION'][:] = cat_lim1['DECAM_MW_TRANSMISSION']
        cat_lim['WISE_MW_TRANSMISSION'][:] = cat['WISE_MW_TRANSMISSION'][mine][mask]

    targets['INTRINSIC_NOISELEVEL'][:, :6] = (cat_lim['DECAM_DEPTH'] ** -0.5 / cat_lim['DECAM_MW_TRANSMISSION'])
    targets['INTRINSIC_NOISELEVEL'][:, 6:] = (cat_lim['WISE_FLUX_IVAR'] ** -0.5 / cat_lim['WISE_MW_TRANSMISSION'])

    nanmask = np.isnan(targets['INTRINSIC_NOISELEVEL'])
    targets['INTRINSIC_NOISELEVEL'][nanmask] = np.inf

    targets['CONFIDENCE'] = targets['INTRINSIC_FLUX'] / targets['INTRINSIC_NOISELEVEL']

    targets = np.concatenate(comm.allgather(targets))

    if comm.rank == 0:
        print('Total number of objects selected', len(targets))

    return targets

if __name__=="__main__":

    targets = select_objs(decals, ns)

    if MPI.COMM_WORLD.rank == 0:
        with h5py.File(ns.output, 'w') as ff:
            ds = ff.create_dataset('_HEADER', shape=(0,))
            for key, value in ns.__dict__.items():
                ds.attrs[key] = value
            for column in targets.dtype.names:
                ds = ff.create_dataset(column, data=targets[column])
