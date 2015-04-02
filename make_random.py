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

from __future__ import print_function,division

__author__ = "Yu Feng and Martin White"
__version__ = "1.0"
__email__  = "yfeng1@berkeley.edu or mjwhite@lbl.gov"

import numpy as N
import sys

#from mpi4py            import MPI
from model.utils       import sharedmem
from model.sfdmap      import SFDMap
from model.datarelease import DataRelease
def fill_random(dr, Nran, seed=999993):
    # Generate uniformly distributed points within the boundary (decimal deg).
    # Everyone generates all of the points, then we specialize to a subset
    #

    # use the same seed to ensure the random points are identical cross ranks

    rng = N.random.RandomState(seed)

    ramin,ramax,dcmin,dcmax = dr.footprint.range

    coord = N.empty((2, Nran))
    start = 0
    while start != Nran:
        with sharedmem.MapReduce() as pool:

            chunksize = 1024 * 512
            u1,u2= rng.uniform(size=(2, 1024 * 1024 * 32))

            def work(i):
                # find my slice
                myu1 = u1[i:i+chunksize]
                myu2 = u2[i:i+chunksize]

                cmin = N.sin(dcmin*N.pi/180)
                cmax = N.sin(dcmax*N.pi/180)

                RA   = ramin + myu1*(ramax-ramin)
                DEC  = 90-N.arccos(cmin+myu2*(cmax-cmin))*180./N.pi

                # filter out those not in any bricks
                coord1 = dr.footprint.filter((RA, DEC))
                
                return coord1

            coord1 = N.concatenate(
                    pool.map(work, range(0, 1024 * 1024 * 32, chunksize)),
                    axis=-1)

        coord1 = coord1[:, :min(len(coord1.T), Nran - start)]
        sl = slice(start, start + len(coord1.T))
        coord[:, sl] = coord1
        start = start + len(coord1.T)
        print(start, '/', Nran, 'filled')

    if len(N.unique(dr.brickindex.query(coord))) != len(dr.footprint.bricks):
        # if this happens, we don't have enough points to cover the foot print
        # fairly. Some bricks have no random points.
        raise RuntimeError("Too few random points, increase Nran")

    return coord

def apply_samp_cut(coord, dr, sfd, samp):
    """ apply the cuts for samp 
        returns ww (indices of good samples), rmag (r band mag)
        
        This assumes every point (assuming LF is a delta function at the faint end cut)
    """
    # Get the flux depths and MW transmission in the relevant filters.
    # Everyone needs "r", but only LRGs need "z" and "W1".
    # For out-of-bounds points, return 0.

    ebv  = sfd.extinction(None, coord[0], coord[1],get_ebv=True)
    rdep = dr.readout(coord,dr.images['depth']['r'],\
                      default=0,ignore_missing=True)
    rtrn = 10.0**(-ebv*dr.extinction['r']/2.5)
    if samp=="LRG":
        zdep = dr.readout(coord,dr.images['depth'][ 'z'],\
                          default=0,ignore_missing=True)
        wdep = dr.readout(coord,dr.images['depth']['W1'],\
                          default=0,ignore_missing=True)
        ztrn = 10.0**(-ebv*dr.extinction[ 'z']/2.5)
        wtrn = 10.0**(-ebv*dr.extinction['W1']/2.5)
    # For now we use a 5-sigma cut in extinction-correct flux as our limit.
    # Recall "depth" is stored as inverse variance.
    rlim = 5.0/N.sqrt(rdep+1e-30) / rtrn
    # It's also useful to have r magnitude later.
    rmag = 22.5-2.5*N.log10( rlim.clip(1e-15,1e15) )
    #
    # Now work out which points pass the cuts--note we want a flux
    # limit *lower* than the catalog cutoff.
    # This code should probably be refactored so that it is
    # connected in some way to the target selection code.
    # For now it just returns a list of indices passing the
    # cuts, so the functionality should be easy to reproduce.
    #
    if samp=="LRG":
        zlim = 5.0/N.sqrt(zdep+1e-30) / ztrn
        wlim = 5.0/N.sqrt(wdep+1e-30) / wtrn
        ww   = N.nonzero( (rlim<10.0**((22.5-23.00)/2.5))&\
                          (zlim<10.0**((22.5-20.56)/2.5))&\
                          (wlim<10.0**((22.5-19.50)/2.5)) )[0]
    elif samp=="ELG":
        ww = N.nonzero( rlim<10.0**((22.5-23.40)/2.5) )[0]
    elif samp=="QSO":
        ww = N.nonzero( rlim<10.0**((22.5-23.00)/2.5) )[0]
    else:
        raise RuntimeError,"Unknown sample "+samp
    return ww, rmag

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

    print ('making randoms in the survey footprint')

    coord = fill_random(dr, Nran)

    print ('querying the depth of randoms')

    # The call to optimize reorders the array to be brick-ordered.
    # read out is faster this way

    coord = dr.brickindex.optimize(coord)

    with sharedmem.MapReduce() as pool:

        chunksize = 1024 * 8

        # purge the file once for all

        fout  = "randoms_%s.rdz"%samp
        with open(fout, 'w') as ff:
            pass
        
        def work(i):
            # Work out which points belong on which slice.
            mycoord = coord[:, i:i+chunksize]

            # apply the cut for sample type on my slice

            ww, rmag = apply_samp_cut(mycoord, dr, sfd, samp)

            # notify user this slide is done
            print(i, '/', len(coord.T))

            return mycoord, ww, rmag

        def reduce(mycoord, ww, rmag):
            #
            # and take turns writing this out:
            #
            with open(fout, 'a') as ff:
                for j in ww:
                    ff.write("%15.10f %15.10f %15.10f %15.10f\n"%\
                      (mycoord[0][j],mycoord[1][j],0.5,rmag[j]))

        pool.map(work, range(0, len(coord.T), chunksize), reduce=reduce)
    print ('done')
    #


if __name__ == '__main__':    
    for samp in ["ELG"]:
        make_random(samp)
    #
