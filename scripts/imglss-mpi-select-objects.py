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
from imaginglss.analysis    import cuts
from imaginglss.model       import dataproduct
from imaginglss.cli         import CLI

cli = CLI("Select Objects based on Target definitions", enable_target_plugins=True)
cli.add_argument("--use-depth-bricks", action='store_true', default=False, help="use tractor's brick depth in the catalogue, very slow!")
cli.add_argument("--use-bigfile", action='store_true', default=False ,help='save as a bigfile; use bigfile-convert to convert afterwards')

cli.add_argument("--limit", type=int, default=None, help='limit to use this many candidates.')
cli.add_target_type_argument("ObjectType")
cli.add_argument("output", help="Output file name. A new object catalogue file will be created.")

ns = cli.parse_args()

decals = DECALS(ns.conf)

from mpi4py import MPI

np.seterr(divide='ignore', invalid='ignore')

def select_objs(decals, ns, comm=MPI.COMM_WORLD):
    """
    Does the actual selection, imposing cuts on the fluxes
    """

    # Get instances of a data release and SFD dust map.
    dr     = decals.datarelease
    sfd    = decals.sfdmap
    cat    = dr.catalogue
    totalsize = cat.size
    if ns.limit is not None:
        totalsize = ns.limit
    #
    mystart = totalsize * comm.rank // comm.size
    myend = totalsize * (comm.rank + 1) // comm.size
    #
    mine = slice(mystart, myend)

    if comm.rank == 0:
        print('Rank 0 with', myend - mystart, 'items', 'total', totalsize)

    with dr.catalogue as cat:
        rows = cat[mine]
        mask = cuts.apply(comm, ns.ObjectType, rows)

    if comm.rank == 0:
        print('Rank 0 selected', mask.sum(), 'items')

    targets = np.empty(mask.sum(), dtype=dataproduct.ObjectCatalogue)

    targets['RA']   = cat[ 'RA'][mine][mask]
    targets['DEC']   = cat['DEC'][mine][mask]

    for i, flux, mw in [
            (1, 'FLUX_G', 'MW_TRANSMISSION_G'),
            (2, 'FLUX_R', 'MW_TRANSMISSION_R'),
            (4, 'FLUX_Z', 'MW_TRANSMISSION_Z'),
            (6, 'FLUX_W1', 'MW_TRANSMISSION_W1'),
            (7, 'FLUX_W2', 'MW_TRANSMISSION_W2')]:

        GFLUX = cat['FLUX_G'][mine][mask] / cat['MW_TRANSMISSION_G'][mine][mask]
        RFLUX = cat['FLUX_R'][mine][mask] / cat['MW_TRANSMISSION_R'][mine][mask]
        ZFLUX = cat['FLUX_Z'][mine][mask] / cat['MW_TRANSMISSION_Z'][mine][mask]
        W1FLUX = cat['FLUX_W1'][mine][mask] / cat['MW_TRANSMISSION_W1'][mine][mask]
        W2FLUX = cat['FLUX_W2'][mine][mask] / cat['MW_TRANSMISSION_W2'][mine][mask]

        targets['INTRINSIC_FLUX'][:, i] = (cat[flux][mine][mask] / cat[mw][mine][mask]).clip(1e-15, 1e15)

    # Now we need to pass this through our mask since galaxies can
    # appear even in regions where our nominal depth is insufficient
    # for a complete ns.ObjectType.

    # note that although DR2 has these numbers, we query the raw images
    # to keep it consistent with make_randoms.py.  This shall agree with
    # the DR2 tractor file numbers.

    # ... and extract only the objects which passed the cuts.

    # we convert cat-lim to the DR3-like schema, because read_depth is doing
    # that.
    cat_lim = np.empty(len(targets), dtype=[
        ('DECAM_DEPTH', ('f4', 6)),
        ('WISE_FLUX_IVAR', ('f4', 4)),
        ('DECAM_MW_TRANSMISSION', ('f4', 6)),
        ('WISE_MW_TRANSMISSION', ('f4', 4)),
        ])

    if not ns.use_depth_bricks:
        for i, depth, mw in [
            (1, 'PSFDEPTH_G', 'MW_TRANSMISSION_G'),
            (2, 'PSFDEPTH_R', 'MW_TRANSMISSION_R'),
            (4, 'PSFDEPTH_Z', 'MW_TRANSMISSION_Z')
        ]:

            cat_lim['DECAM_DEPTH'][:, i] = cat[depth][mine][mask]
            cat_lim['DECAM_MW_TRANSMISSION'][:, i] = cat[mw][mine][mask]
    else:
        cat_lim1 = dr.read_depths((targets['RA'], targets['DEC']), 'grz')
        cat_lim['DECAM_DEPTH'][:, i] = cat_lim1['DECAM_DEPTH']
        cat_lim['DECAM_MW_TRANSMISSION'][:, i] = cat_lim1['DECAM_MW_TRANSMISSION']

    for i, ivar, mw in [
        (0, 'FLUX_IVAR_W1', 'MW_TRANSMISSION_W1'),
        (1, 'FLUX_IVAR_W2', 'MW_TRANSMISSION_W2')
        ]:
        cat_lim['WISE_FLUX_IVAR'][:, i] = cat[ivar][mine][mask]
        cat_lim['WISE_MW_TRANSMISSION'][:, i] = cat[mw][mine][mask]

    targets['INTRINSIC_NOISELEVEL'][:, :6] = (cat_lim['DECAM_DEPTH'] ** -0.5 / cat_lim['DECAM_MW_TRANSMISSION'])
    targets['INTRINSIC_NOISELEVEL'][:, 6:] = (cat_lim['WISE_FLUX_IVAR'] ** -0.5 / cat_lim['WISE_MW_TRANSMISSION'])

    nanmask = np.isnan(targets['INTRINSIC_NOISELEVEL'])
    targets['INTRINSIC_NOISELEVEL'][nanmask] = np.inf

    targets['CONFIDENCE'] = targets['INTRINSIC_FLUX'] / targets['INTRINSIC_NOISELEVEL']

    return targets

def savehdf(ns, targets, comm):
    size = comm.allreduce(len(targets))

    if comm.rank == 0:
        print('Total number of objects selected', size)

        with h5py.File(ns.output, 'w') as ff:
            ds = ff.create_dataset('_HEADER', shape=(0,))
            ds.attrs.update(cli.prune_namespace(ns))

            for column in targets.dtype.names:
                shape = tuple([size] + list(targets[column].shape[1:]))
                dtype = targets[column].dtype
                ds = ff.create_dataset(column, shape=shape, dtype=dtype)

    comm.barrier()

    # the loop makes sure the ranks take turns; avoiding race condition
    # with HDF5. Really could have saved as a bigfile instead to avoid
    # this craziness.

    offset = sum(comm.allgather(len(targets))[:comm.rank])

    for i in range(comm.size):
        comm.barrier()

        if i != comm.rank :
            continue

        print('rank', comm.rank, 'writing')
        for column in targets.dtype.names:
            print('rank', comm.rank, 'writing', column)
            with h5py.File(ns.output, 'r+') as ff:
                print('rank', comm.rank, 'writing', column, ff)
                for column in targets.dtype.names:
                    ff[column][offset:offset+len(targets)] = targets[column]

def savebigfile(ns, targets, comm):
    import bigfile
    with bigfile.BigFileMPI(comm, ns.output, create=True) as ff:
        for column in targets.dtype.names:
            with ff.create_from_array(column, targets[column], Nfile=1) as bb:
                pass

if __name__=="__main__":

    targets = select_objs(decals, ns)

    comm = MPI.COMM_WORLD
    if not ns.use_bigfile:
        savehdf(ns, targets, comm)
    else:
        savebigfile(ns, targets, comm)

