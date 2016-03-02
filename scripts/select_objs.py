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
from imaginglss             import DECALS
from imaginglss.analysis    import targetselection
from imaginglss.analysis    import completeness
from imaginglss.analysis    import tycho_veto
from imaginglss.analysis    import cuts
from imaginglss.utils       import output

from argparse import ArgumentParser

ap = ArgumentParser("select_objs.py")
ap.add_argument("ObjectType", choices=[i for i in targetselection.__all__])
ap.add_argument("output", type=output.writer)
ap.add_argument("--use-tractor-depth", action='store_true', default=False, help="Use Tractor's Depth in the catalogue")
ap.add_argument("--sigma-z", type=float, default=3.0)
ap.add_argument("--sigma-g", type=float, default=5.0)
ap.add_argument("--sigma-r", type=float, default=5.0)
ap.add_argument("--with-tycho", choices=[i for i in dir(tycho_veto) if not str(i).startswith( '_' )], help="Type of veto.")
ap.add_argument("--conf", default=None,
        help="Path to the imaginglss config file, default is from DECALS_PY_CONFIG")

ns = ap.parse_args()
ns.conf = DECALS(ns.conf)

from mpi4py import MPI

np.seterr(divide='ignore', invalid='ignore')

def select_objs(ns, comm=MPI.COMM_WORLD):
    """
    Does the actual selection, imposing cuts on the fluxes
    """

    dtype = np.dtype([
        ('RA', 'f8'),
        ('DEC', 'f8'),
        ('PHOTO_Z', 'f8'),
        ('DECAM_MAG', ('f4', 6)),
        ('DECAM_NOISE_LEVEL', ('f4', 6)),
        ('COMPLETE','f8'),
        ])

    # Get instances of a data release and SFD dust map.
    decals = ns.conf
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

    CANDIDATES = np.empty(mask.sum(), dtype=dtype)
    CANDIDATES['PHOTO_Z']   = -1.0
    CANDIDATES['RA']   = cat[ 'RA'][mine][mask]
    CANDIDATES['DEC']   = cat['DEC'][mine][mask]
 
    FLUX = cat['DECAM_FLUX'][mine][mask] / cat['DECAM_MW_TRANSMISSION'][mine][mask]
    CANDIDATES['DECAM_MAG'] = 22.5-2.5*np.log10((FLUX).clip(1e-15,1e15))

    # Now we need to pass this through our mask since galaxies can
    # appear even in regions where our nominal depth is insufficient
    # for a complete ns.ObjectTypee.

    # note that although DR2 has these numbers, we query the raw images
    # to keep it consistent with make_randoms.py.  This shall agree with
    # the DR2 tractor file numbers.

    # ... and extract only the objects which passed the cuts.

    compcut = getattr(completeness, ns.ObjectType)
    sigma = dict(z=ns.sigma_z, g=ns.sigma_g, r=ns.sigma_r)
    if ns.use_tractor_depth:
        cat_lim = np.empty(len(CANDIDATES), dtype=[
            ('DECAM_DEPTH', cat.dtype['DECAM_DEPTH']),
            ('DECAM_MW_TRANSMISSION', cat.dtype['DECAM_MW_TRANSMISSION']),
            ])
        cat_lim['DECAM_DEPTH'][:] = cat['DECAM_DEPTH'][mine][mask]
        cat_lim['DECAM_MW_TRANSMISSION'][:] = cat['DECAM_MW_TRANSMISSION'][mine][mask]
    else:
        cat_lim = dr.read_depths((CANDIDATES['RA'], CANDIDATES['DEC']), compcut.bands)

    # It's also useful to store the 1 sigma limits for later
    CANDIDATES['DECAM_NOISE_LEVEL'] = (cat_lim['DECAM_DEPTH'] ** -0.5 / cat_lim['DECAM_MW_TRANSMISSION']).clip(0, 60)

    for band in compcut.bands:
        ind = dr.bands[band]
        missing_depth = sum(comm.allgather(
                (cat_lim['DECAM_DEPTH'][:, ind] == 0).sum()))
        if comm.rank == 0:
            print('Objects in bricks with missing depth images (',band,'): ',\
                  missing_depth)

    mask = cuts.apply(comm, compcut(sigma), cat_lim)

    total_complete = sum(comm.allgather(mask.sum()))
    if comm.rank == 0:
        print('Total number of objects in complete area: ', total_complete)

    
    CANDIDATES['COMPLETE'] = 0
    CANDIDATES['COMPLETE'][mask] = 1
    
    if ns.with_tycho is not None:
        veto = getattr(tycho_veto, ns.with_tycho)
        mask = veto(decals.tycho, (CANDIDATES['RA'], CANDIDATES['DEC']))

        total_complete = sum(comm.allgather(mask.sum()))
        if comm.rank == 0:
            print('Using tycho veto', ns.with_tycho,'...')
            print('Total number of objects not close to stars', total_complete)

        CANDIDATES = CANDIDATES[mask]
    else:
        if comm.rank == 0:
            print('Not applying cuts for star proximity.')
        

    CANDIDATES  = np.concatenate(comm.allgather(CANDIDATES))

    return CANDIDATES

if __name__=="__main__":

    CANDIDATES = select_objs(ns)

    if MPI.COMM_WORLD.rank == 0:
        ns.output.write(CANDIDATES, ns.__dict__)
