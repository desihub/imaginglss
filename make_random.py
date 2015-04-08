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
# Alternatively we could impose the cuts before writing the
# files, thus generating specific randoms "at once".
# This code takes the latter approach, while "mpirandom" takes
# the former approach.
#
# Currently this code uses OpenMP parallelization, the switch to
# using MPI is relatively straightforward.
#
from __future__ import print_function,division

__author__ = "Yu Feng and Martin White"
__version__ = "1.0"
__email__  = "yfeng1@berkeley.edu or mjwhite@lbl.gov"

import numpy             as N
import select_elgs       as S
from   model.utils       import sharedmem
from   model.sfdmap      import SFDMap
from   model.datarelease import DataRelease
import sys
#from mpi4py            import MPI

import cuts


def fill_random(dr,Nran,seed=999993):
    """
    fill_random(dr,Nran,seed=999993):
    Generate uniformly distributed points within the boundary that lie in
    bricks.  We generate in the ra/dec area, then remove points not in any
    bricks.  This hugely increases the memory efficiency for footprints,
    like DR1, where the survey covers several disjoint patches of the sky.
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
        print(start, '/', Nran, 'filled')
    if len(N.unique(dr.brickindex.query(coord))) != len(dr.footprint.bricks):
        # If this happens, we don't have enough points to cover the foot print
        # fairly. Some bricks have no random points.
        raise RuntimeError("Too few random points, increase Nran")
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
    # Get the flux depths and MW transmission in the relevant filters.
    # For out-of-bounds points, return 0.

    ebv  = sfd.ebv(coord[0], coord[1])

    sigma = 5
    # Now work out which points pass the cuts -- note we want a flux
    # limit *lower* than the catalog cutoff.
    # This code should be related to the target selection code, for now
    # we are hard-wiring it and calling the similar routine in select_elgs.
    # The code just returns a list of indices passing the cuts, so the
    # functionality should be easy to reproduce.
    if samp=="LRG":
        rlim, zlim, wlim = S.findlim(dr, sfd, coord, ['r', 'z', 'W1'], sigma=5)

        mask = cuts.Completeness.LRG(rlim, zlim, wlim)
    elif samp=="ELG":
        rlim, glim, zlim = S.findlim(dr, sfd, coord, ['r', 'g', 'z'], sigma=5)
        mask = cuts.Completeness.ELG(rlim, glim, zlim)

    elif samp=="QSO":
        rlim, = S.findlim(dr, sfd, coord, ['r'], sigma=5)
        mask = cuts.Completeness.QSO(rlim)
    else:
        raise RuntimeError,"Unknown sample "+samp
    ww = N.nonzero(mask)[0]
    # It's also useful to have r magnitude later.
    rmag = 22.5-2.5*N.log10( rlim.clip(1e-15,1e15) )
    return( (ww,rmag) )
    #



def make_random(samp,Nran=10000000):
    """
    make_random(samp,Nran=1000000)
    Does the work of making randoms.  The sample type is passed as a string.
    """
    if samp not in ["LRG","ELG","QSO"]:
        raise RuntimeError,"Unknown sample "+samp
    # Get the total footprint bounds, to throw randoms within, and an E(B-V)
    # map instance.
    dr = DataRelease(version='DR1') #, _covered_brickids=subset)
    sfd= SFDMap(dustdir="/project/projectdirs/desi/software/edison/dust/v0_0/")
    print('Making randoms in the survey footprint.')
    coord = fill_random(dr, Nran)
    print('Querying the depth of randoms')

    # The call to optimize reorders the array to be brick-ordered.
    # Read out is faster this way
    coord = dr.brickindex.optimize(coord)

    bid = dr.brickindex.query(coord)
    print('unique bricks', len(N.unique(bid)))

    with sharedmem.MapReduce() as pool:
        # prepare a parallel section

        chunksize = max(len(coord.T) // (pool.np * 4), 1)

        # purge the file once for all
        fout  = "randoms_%s.rdz"%samp
        with open(fout,'w') as ff:
            pass
        def work(i):
            print (i)
            # Work out which points belong on which slice.
            mycoord = coord[:, i:i+chunksize]
            # Apply the cut for sample type on my slice
            ww,rmag = apply_samp_cut(mycoord,dr,sfd,samp)
            # Notify user this slice is done
            print(i, '/', len(coord.T))
            return( (mycoord,ww,rmag) )
        def reduce(mycoord, ww, rmag):
            # and take turns writing this out:
            with open(fout, 'a') as ff:
                for j in ww:
                    ff.write("%15.10f %15.10f %15.10f %15.10f\n"%\
                      (mycoord[0][j],mycoord[1][j],0.5,rmag[j]))
            return len(ww)
        total = sum(pool.map(work,range(0,len(coord.T),chunksize),reduce=reduce))
    print('Total area (sq.deg.) ',dr.footprint.area*1.0*total/len(coord.T))
    print('Done!')
    #


 
if __name__ == '__main__':    
    for samp in ["ELG"]:
        make_random(samp, Nran=1000000)
    #
