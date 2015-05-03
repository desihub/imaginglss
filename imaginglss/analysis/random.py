import numpy as N
from ..utils import sharedmem

def fill(dr,Nran,seed=999993, verbose=False):
    """
    Generate uniformly distributed points within the boundary that lie in
    bricks.  

    We generate in the ra/dec area, then remove points not in any
    bricks.  This hugely increases the memory efficiency for footprints,
    like DR1, where the survey covers several disjoint patches of the sky.

    Parameters
    ----------
    dr  : :py:class:`model.datarelease.DataRelease`
        the data release objects
    Nran  : int
        number of random points to generate
    seed : int
        random seed
    verbose : boolean
        display progress to stdout. Do not use.

    Returns
    -------
    coord   : array_like (2, N)
        coord = (RA, DEC) is the coordinate of the random positions with in the
        footprint of dr.

    Raises
    ------
    RuntimeError :
        if the number of random points is so few it won't cover the footprint
        with at least 1 point per brick.

    """
    rng = N.random.RandomState(seed)
    ramin,ramax,dcmin,dcmax = dr.footprint.range
    coord = N.empty((2, Nran))
    start = 0
    while start != Nran:
        with sharedmem.MapReduce() as pool:
            # prepare for a parallel section.
            chunksize = 1024 * 512
            u1,u2= rng.uniform(size=(2, 1024 * 1024 * 32))
            def work(i):
                # Find my slice
                myu1 = u1[i:i+chunksize]
                myu2 = u2[i:i+chunksize]
                #
                cmin = N.sin(dcmin*N.pi/180)
                cmax = N.sin(dcmax*N.pi/180)
                #
                RA   = ramin + myu1*(ramax-ramin)
                DEC  = 90-N.arccos(cmin+myu2*(cmax-cmin))*180./N.pi
                # Filter out those not in any bricks: only very few points remain
                coord1 = dr.footprint.filter((RA, DEC))
                return(coord1)
            # Run; the number here is just make sure each batch won't 
            # use too much memory
            coord1 = N.concatenate(
                    pool.map(work, range(0,1024*1024*32,chunksize)),
                    axis=-1)
        # Are we full?
        coord1 = coord1[:, :min(len(coord1.T), Nran - start)]
        sl = slice(start, start + len(coord1.T))
        coord[:, sl] = coord1
        start = start + len(coord1.T)
        if verbose:
            print(start, '/', Nran, 'filled')
    if len(N.unique(dr.brickindex.query(coord))) != len(dr.footprint.bricks):
        # If this happens, we don't have enough points to cover the foot print
        # fairly. Some bricks have no random points.
        raise RuntimeError("Too few random points, increase Nran")
    return coord
