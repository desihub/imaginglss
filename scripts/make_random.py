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

from imaginglss.analysis    import tycho_veto

from argparse import ArgumentParser

ap = ArgumentParser("make_random.py")
ap.add_argument("Nran", type=int, help="Minimum number of randoms")
ap.add_argument("ObjectType", choices=["QSO", "LRG", "ELG", "BGS"])
ap.add_argument("output")
ap.add_argument("--sigma-z", type=float, default=3.0)
ap.add_argument("--sigma-g", type=float, default=5.0)
ap.add_argument("--sigma-r", type=float, default=5.0)
ap.add_argument("--with-tycho", choices=[i for i in dir(tycho_veto) if not str(i).startswith( '_' )], help="Type of veto.")
ap.add_argument("--conf", default=None, 
        help="Path to the imaginglss config file, default is from DECALS_PY_CONFIG")

ns = ap.parse_args()

import numpy             as np
from   imaginglss             import DECALS
from   imaginglss.model.datarelease import Footprint
from mpi4py            import MPI

from   imaginglss.analysis    import cuts
from imaginglss.analysis    import completeness
from imaginglss.analysis    import tycho_veto

np.seterr(divide='ignore', invalid='ignore')

def fill_random(footprint, Nran, rng):
    """
    Generate uniformly distributed points within the boundary that lie in
    bricks.  We generate in the ra/dec area, then remove points not in any
    bricks.  This hugely increases the memory efficiency for footprints,
    like DR1, where the survey covers several disjoint patches of the sky.

    """


    coord = np.empty((2, Nran))

    ramin,ramax,dcmin,dcmax,area = footprint.range

    start = 0
    while start != Nran:
        # prepare for a parallel section.
        chunksize = 1024 * 512
        u1,u2= rng.uniform(size=(2, 1024 * 1024) )

        #
        cmin = np.sin(dcmin*np.pi/180)
        cmax = np.sin(dcmax*np.pi/180)
        #
        RA   = ramin + u1*(ramax-ramin)
        DEC  = 90-np.arccos(cmin+u2*(cmax-cmin))*180./np.pi
        # Filter out those not in any bricks: only very few points remain
        coord1 = footprint.filter((RA, DEC))

        # Are we full?
        coord1 = coord1[:, :min(len(coord1.T), Nran - start)]
        sl = slice(start, start + len(coord1.T))
        coord[:, sl] = coord1
        start = start + len(coord1.T)
    
    return coord


def make_random(ns, comm=MPI.COMM_WORLD):
    """
    Does the work of making randoms.  The ns.ObjectTypee type is passed as a string.
    """

    # Get the total footprint bounds, to throw randoms within, and an E(B-V)
    # map instance.
    decals = DECALS(ns.conf)
    dr = decals.datarelease
    sfd= decals.sfdmap

    rng = np.random.RandomState(99934123)
    seeds = rng.randint(99999999, size=comm.size)
    rng = np.random.RandomState(seeds[comm.rank])
    if comm.rank == 0:
        print('Making randoms in the survey footprint.')

    # create a smaller footprint for just this rank

    Nbricks = len(dr.footprint.bricks)
    mybricks = dr.footprint.bricks[comm.rank * Nbricks // comm.size
                     : (comm.rank + 1) * Nbricks // comm.size]
    footprint = Footprint(mybricks, dr.brickindex)

    Nbricks = comm.allgather(len(mybricks))

    # Distribute myns.Nran, proportional to the number of bricks
    # important if a rank has no bricks, it expects no randoms

    myNran = rng.poisson(footprint.area / dr.footprint.area * ns.Nran)
    
    print (comm.rank, 'has', len(mybricks), 'bricks', myNran, 'randoms')
    # fill it with random points 
    coord = fill_random(footprint, myNran, rng)

    ns.Nran = sum(comm.allgather(len(coord[0])))

    #Currently this applies cuts for the "100%" complete ns.ObjectTypee, i.e. assuming
    #the LF is a delta function at the faint end cut.  This avoids needing
    #to know the LF for now.  Later we could keep a fraction of the randoms
    #based on the survey limit (and the LF) or we could produce randoms for
    #"100%" complete ns.ObjectTypees of different luminosity thresholds.

    compcut = getattr(completeness, ns.ObjectType)
    sigma = dict(z=ns.sigma_z, g=ns.sigma_g, r=ns.sigma_r)
    cat_lim = dr.read_depths(coord, compcut.bands)

    mask = cuts.apply(comm, compcut(sigma), cat_lim)

    # It's also useful to have r magnitude later.
    LIM = cat_lim['DECAM_FLUX_IVAR'] ** -0.5 / cat_lim['DECAM_MW_TRANSMISSION']

    total_complete = sum(comm.allgather(mask.sum()))
    if comm.rank == 0:
        print('Total number of objects in complete area: ', total_complete)
    RA, DEC = coord
    RA  = RA [mask]
    DEC = DEC[mask]
    LIM = LIM[mask]

    if ns.with_tycho is not None:
        veto = getattr(tycho_veto, ns.with_tycho)
        mask = veto(decals.tycho, (RA, DEC))

        total_complete = sum(comm.allgather(mask.sum()))
        if comm.rank == 0:
            print('Total number of objects not close to stars', total_complete)

        RA  = RA [mask]
        DEC = DEC[mask]
        LIM = LIM[mask]
    else:
        if comm.rank == 0:
            print('Not applying cuts for star proximity.')

    LIM = comm.gather(LIM)
    coord = comm.gather((RA, DEC))
    if comm.rank == 0:
        coord = np.concatenate(coord, axis=-1)
        LIM = np.concatenate(LIM)
        fraction = len(coord[0]) * 1.0 / ns.Nran
        print('Accept rate', fraction)
        print('Total area (sq.deg.) ',dr.footprint.area * fraction)
        print('Done!')
    return coord, LIM
    #


if __name__ == '__main__':    

    coord, LIM = make_random(ns)

    if MPI.COMM_WORLD.rank == 0:
        with open(ns.output,'w') as ff:
            ff.write("# sigma_z=%g sigma_g=%g sigma_r=%g\n" % (ns.sigma_z, ns.sigma_g, ns.sigma_r))
            ff.write("# ra dec weight ")
            for band in 'ugrizY':
                ff.write(" %15s" % band)
            ff.write("\n")
            
            for j in range(0, len(LIM)):
                ff.write("%15.10f %15.10f %15.10f "%\
                  (coord[0][j],coord[1][j],0.5))
                for band in range(6):
                    ff.write(" %15.10f" % LIM[j][band])
                ff.write("\n")

