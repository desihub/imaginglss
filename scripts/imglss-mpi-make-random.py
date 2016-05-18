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

import h5py

from   imaginglss             import DECALS

from imaginglss.cli import CLI

cli = CLI("Generate Uniform Randoms in the footprint and query the depth.")
cli.add_argument("Nran", type=int, help="Minimum number of randoms")
cli.add_argument("output", help="File to store the randoms. Will be created." )

ns = cli.parse_args()
decals = DECALS(ns.conf)

import numpy             as np
from   imaginglss.model.datarelease import Footprint
from   imaginglss.model import dataproduct
from   mpi4py            import MPI
from   imaginglss.analysis    import cuts

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


def make_random(decals, ns, comm=MPI.COMM_WORLD):
    """
    Does the work of making randoms.  The ns.ObjectTypee type is passed as a string.
    """

    # Get the total footprint bounds, to throw randoms within, and an E(B-V)
    # map instance.
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

    if comm.rank == 0:
        print('Not applying cuts for star proximity.')

    randoms = comm.gather(randoms)
    if comm.rank == 0:
        randoms = np.concatenate(randoms)

        print('Total number of objects ', len(randoms))
        fraction = len(randoms) * 1.0 / ns.Nran
        print('Accept rate', fraction)
        print('Total area (sq.deg.) ',dr.footprint.area * fraction)
        print('Done!')
    return randoms
    
if __name__ == '__main__':    

    randoms = make_random(decals, ns)

    if MPI.COMM_WORLD.rank == 0:
        with h5py.File(ns.output, 'w') as ff:
            ds = ff.create_dataset('_HEADER', shape=(0,))
            for key, value in ns.__dict__.items():
                ds.attrs[key] = value
            for column in randoms.dtype.names:
                ds = ff.create_dataset(column, data=randoms[column])
