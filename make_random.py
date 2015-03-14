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

from mpi4py            import MPI
from model.sfdmap      import SFDMap
from model.datarelease import DataRelease


def make_random(samp,Nran=1000000):
    """
    make_random(samp,Nran=1000000):
    Does the work of making randoms.  The sample type is passed as a string.
    """
    if samp not in ["LRG","ELG","QSO"]:
        raise RuntimeError,"Unknown sample "+samp
    # Get the total footprint bounds, to throw randoms within, and an E(B-V)
    # map instance.
    dr = DataRelease(version='EDR')
    ramin,ramax,dcmin,dcmax = dr.observed_range
    sfd= SFDMap(dustdir="/project/projectdirs/desi/software/edison/dust/v0_0/")
    # Generate uniformly distributed points within the boundary (decimal deg).
    # Everyone generates all of the points, then we specialize to a subset.
    # The call to optimize reorders the array to be brick-ordered.
    u1,u2= N.random.uniform(size=(2,Nran))
    cmin = N.sin(dcmin*N.pi/180)
    cmax = N.sin(dcmax*N.pi/180)
    RA   = ramin + u1*(ramax-ramin)
    DEC  = 90-N.arccos(cmin+u2*(cmax-cmin))*180./N.pi
    (RA,DEC) = dr.brickindex.optimize((RA,DEC))
    # Work out which points belong on which slice.
    comm    = MPI.COMM_WORLD
    mystart =  comm.rank   *len(RA)//comm.size
    myend   = (comm.rank+1)*len(RA)//comm.size
    sl      = slice(mystart,myend)
    coord   = (RA[sl],DEC[sl])
    # Get the flux depths and MW transmission in the relevant filters.
    # Everyone needs "r", but only LRGs need "z" and "W1".
    # For out-of-bounds points, return 0.
    ebv  = sfd.extinction(None,RA[sl],DEC[sl],get_ebv=True)
    rdep = dr.readout(coord,dr.images['depth']['r'],default=0)
    rtrn = 10.0**(-ebv*dr.extinction['r']/2.5)
    if samp=="LRG":
        zdep = dr.readout(coord,dr.images['depth'][ 'z'],default=0)
        wdep = dr.readout(coord,dr.images['depth']['W1'],default=0)
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
    #
    # and take turns writing this out:
    fout  = "randoms_%s.rdz"%samp
    for i in range(comm.size):
        if i==comm.rank:
            if i==0:
                ff = open(fout,"w")
            else:
                ff = open(fout,"a")
            for j in ww:
                ff.write("%15.10f %15.10f %15.10f %15.10f\n"%\
                  (coord[0][j],coord[1][j],0.5,rmag[j]))
            ff.close()
        comm.Barrier()
    #


if __name__ == '__main__':    
    for samp in ["ELG"]:
        make_random(samp)
    #
