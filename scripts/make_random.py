#!/usr/bin/env python
#
# Generates random points within the footprint based
# on a chosen selection (e.g. LRG, ELG or QSO).
#
# There are two ways we could do this ... first we could
# generate a large number of randoms and store the
# properties required for target selection per random.
# Then a query on this file would return randoms appropriate
# for whatever selection is desired.
#
from __future__ import print_function,division

__author__ = "Yu Feng and Martin White"
__version__ = "1.0"
__email__  = "yfeng1@berkeley.edu or mjwhite@lbl.gov"

import os.path; import sys; sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..'))

from argparse import ArgumentParser

ap = ArgumentParser("make_random.py")
ap.add_argument("Nran", type=int, help="Minimum number of randoms")
ap.add_argument("ObjectType", choices=["QSO", "LRG", "ELG", "BGS"])
ap.add_argument("output")
ap.add_argument("--conf", default=None, 
        help="Path to the imaginglss config file, default is from DECALS_PY_CONFIG")

ns = ap.parse_args()

import numpy             as N
from   imaginglss             import DECALS
from   imaginglss.analysis    import cuts
from   imaginglss.model.datarelease import Footprint
from mpi4py            import MPI

N.seterr(divide='ignore', invalid='ignore')

def fill_random(footprint, Nran, rng):
    """
    Generate uniformly distributed points within the boundary that lie in
    bricks.  We generate in the ra/dec area, then remove points not in any
    bricks.  This hugely increases the memory efficiency for footprints,
    like DR1, where the survey covers several disjoint patches of the sky.

    """


    coord = N.empty((2, Nran))

    ramin,ramax,dcmin,dcmax = footprint.range

    start = 0
    while start != Nran:
        # prepare for a parallel section.
        chunksize = 1024 * 512
        u1,u2= rng.uniform(size=(2, 1024 * 1024) )

        #
        cmin = N.sin(dcmin*N.pi/180)
        cmax = N.sin(dcmax*N.pi/180)
        #
        RA   = ramin + u1*(ramax-ramin)
        DEC  = 90-N.arccos(cmin+u2*(cmax-cmin))*180./N.pi
        # Filter out those not in any bricks: only very few points remain
        coord1 = footprint.filter((RA, DEC))

        # Are we full?
        coord1 = coord1[:, :min(len(coord1.T), Nran - start)]
        sl = slice(start, start + len(coord1.T))
        coord[:, sl] = coord1
        start = start + len(coord1.T)
    
    return coord


def apply_samp_cut(coord, dr, sfd, samp):
    """
    apply_samp_cut(coord, dr, sfd, samp): 
    Apply the cuts for sample "samp".  This logic is sufficiently complex it
    is worth having in its own routine.
    Returns ww (indices of good samples), rmag (r band mag)
    Currently this applies cuts for the "100%" complete sample, i.e. assuming
    the LF is a delta function at the faint end cut.  This avoids needing
    to know the LF for now.  Later we could keep a fraction of the randoms
    based on the survey limit (and the LF) or we could produce randoms for
    "100%" complete samples of different luminosity thresholds.
    """

    cut = getattr(cuts.Completeness, samp)

    lim = cuts.findlim(dr,sfd,coord, cut.bands)

    mask = cut(**lim)

    # It's also useful to have r magnitude later.
    rmag = 22.5-2.5*N.log10( lim['r'].clip(1e-15,1e15) )

    return mask, rmag, cut


def make_random(samp, Nran, configfile, output, comm=MPI.COMM_WORLD):
    """
    Does the work of making randoms.  The sample type is passed as a string.
    """

    # Get the total footprint bounds, to throw randoms within, and an E(B-V)
    # map instance.
    decals = DECALS()
    dr = decals.datarelease
    sfd= decals.sfdmap

    rng = N.random.RandomState(99934123)
    seeds = rng.randint(99999999, size=comm.size)
    rng = N.random.RandomState(seeds[comm.rank])
    if comm.rank == 0:
        print('Making randoms in the survey footprint.')

    # create a smaller footprint for just this rank

    Nbricks = len(dr.footprint.bricks)
    mybricks = dr.footprint.bricks[comm.rank * Nbricks // comm.size
                     : (comm.rank + 1) * Nbricks // comm.size]
    footprint = Footprint(mybricks, dr.brickindex)

    Nbricks = comm.allgather(len(mybricks))

    # Distribute myNran, proportional to the number of bricks
    # important if a rank has no bricks, it expects no randoms

    myNran = rng.poisson(footprint.area / dr.footprint.area * Nran)
    
    print (comm.rank, 'has', len(mybricks), 'bricks', myNran, 'randoms')
    # fill it with random points 
    coord = fill_random(footprint, myNran, rng)

    Nran = sum(comm.allgather(len(coord[0])))

    mask, rmag, cut = apply_samp_cut(coord, dr, sfd, samp)

    selected_fraction = 1.0 * sum(comm.allgather(mask.sum(axis=-1))) / Nran 

    if comm.rank == 0:
        print("Selected fraction:\n",
              "\n".join([
                "%s : %g" % (c, f)
                for c, f in zip(cut, selected_fraction)
                ])
            )

    mask = mask.all(axis=0)
    coord = coord[:, mask]
    rmag = rmag[mask]

    coord = comm.gather(coord)
    rmag = comm.gather(rmag)

    if comm.rank == 0:
        coord = N.concatenate(coord, axis=-1)
        rmag = N.concatenate(rmag)

        with open(output,'w') as ff:
            ff.write("# ra dec weight rmag\n")
            for j in range(0, len(rmag)):
                ff.write("%15.10f %15.10f %15.10f %15.10f\n"%\
                  (coord[0][j],coord[1][j],0.5,rmag[j]))

        fraction = len(coord[0]) * 1.0 / Nran
        print('Accept rate', fraction)
        print('Total area (sq.deg.) ',dr.footprint.area * fraction)
        print('Done!')
    #


if __name__ == '__main__':    
    make_random(ns.ObjectType, configfile=ns.conf, Nran=ns.Nran, output=ns.output)
    #
