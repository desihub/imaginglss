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

import os.path; import sys; sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..'))

import numpy as np
from imaginglss             import DECALS
from imaginglss.analysis    import targetselection
from imaginglss.analysis    import completeness
from imaginglss.analysis    import tycho_veto
from imaginglss.analysis    import cuts

from mpi4py import MPI

from argparse import ArgumentParser

ap = ArgumentParser("select_objs.py")
ap.add_argument("ObjectType", choices=["QSO", "LRG", "ELG", "BGS"])
ap.add_argument("output")
ap.add_argument("--sigma-z", type=float, default=3.0)
ap.add_argument("--sigma-g", type=float, default=5.0)
ap.add_argument("--sigma-r", type=float, default=5.0)
ap.add_argument("--with-tycho", help='path to the tycho.fit file for applying veto around bright stars')
ap.add_argument("--conf", default=None, 
        help="Path to the imaginglss config file, default is from DECALS_PY_CONFIG")

ns = ap.parse_args()

np.seterr(divide='ignore', invalid='ignore')

def select_objs(ns, comm=MPI.COMM_WORLD):
    """
    Does the actual selection, imposing cuts on the fluxes
    """
    # Get instances of a data release and SFD dust map.
    decals = DECALS(ns.conf)
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

    # ... and extract only the objects which passed the cuts.
    RA   = cat[ 'RA'][mine][mask]
    DEC   = cat['DEC'][mine][mask]
    FLUX = cat['DECAM_FLUX'][mine][mask] / cat['DECAM_MW_TRANSMISSION'][mine][mask]
    MAG = 22.5-2.5*np.log10((FLUX).clip(1e-15,1e15))

    # Now we need to pass this through our mask since galaxies can
    # appear even in regions where our nominal depth is insufficient
    # for a complete ns.ObjectTypee.
    compcut = getattr(completeness, ns.ObjectType)
    sigma = dict(z=ns.sigma_z, g=ns.sigma_g, r=ns.sigma_r)
    cat_lim = dr.read_depths((RA, DEC), compcut.bands)

    for band in compcut.bands:
        ind = dr.bands[band]
        missing_depth = sum(comm.allgather(
                (cat_lim['DECAM_INVVAR'][:, ind] == 0).sum()))
        if comm.rank == 0:
            print('Objects in bricks with missing depth images (',band,'): ',\
                  missing_depth)

    mask = cuts.apply(comm, compcut(sigma), cat_lim)

    total_complete = sum(comm.allgather(mask.sum()))
    if comm.rank == 0:
        print('Total number of objects in complete area: ', total_complete)

    RA  = RA [mask]
    DEC = DEC[mask]
    MAG = MAG[mask]
    
    if ns.with_tycho is not None:
        tycho = tycho_veto.Tycho(ns.with_tycho)
        veto = getattr(tycho_veto, 'DECAM_' + ns.ObjectType)
        mask = veto(tycho, (RA, DEC))

        total_complete = sum(comm.allgather(mask.sum()))
        if comm.rank == 0:
            print('Total number of objects not close to stars', total_complete)

        RA  = RA [mask]
        DEC = DEC[mask]
        MAG = MAG[mask]
    else:
        if comm.rank == 0:
            print('Not applying cuts for star proximity.')
        

    RA  = np.concatenate(comm.allgather(RA ))
    DEC = np.concatenate(comm.allgather(DEC))
    MAG = np.concatenate(comm.allgather(MAG), axis=0)

    return( (RA,DEC,MAG) )
    #

if __name__=="__main__":

    ra,dc,mag = select_objs(ns)

    if MPI.COMM_WORLD.rank == 0:
        # Just write the ns.ObjectTypee to an ascii text file.
        ff = file(ns.output,"w")
        with ff:
            ff.write("# sigma_z=%g sigma_g=%g sigma_r=%g\n" % (ns.sigma_z, ns.sigma_g, ns.sigma_r))
            ff.write("# %13s %15s %15s" % ("RA","DEC","PhotoZ"))
            for band in 'ugrizY':
                ff.write(" %15s" % band)
            ff.write("\n")
            for i in range(ra.size):
                ff.write("%15.10f %15.10f %15.10f"
                    % (ra[i],dc[i],0.5))
                for band in range(6):
                    ff.write(" %15.10f" % mag[i][band])
                ff.write("\n")
        #
