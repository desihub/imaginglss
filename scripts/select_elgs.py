#!/usr/bin/env python
#
# Code to select ELG targets from Tractor catalogs.
# The target selection criteria are described at:
#   https://desi.lbl.gov/trac/wiki/TargetSelection
#
# These catalogs also need to be vetod using the
# depths and masks which are part of Tractor, since
# catalog entries can appear in places where the
# nominal depth is insufficient to be complete.
#
# Usage: python select_elg.py [--plot]
#
from __future__ import print_function


__author__ = "Yu Feng and Martin White"
__version__ = "1.0"
__email__  = "yfeng1@berkeley.edu or mjwhite@lbl.gov"

import os.path; import sys; sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..'))

import numpy as N
from imaginglss             import DECALS
from imaginglss.analysis    import cuts

from mpi4py import MPI


N.seterr(divide='ignore', invalid='ignore')

def select_elgs(comm=MPI.COMM_WORLD):
    """
    Does the actual selection, imposing cuts on the fluxes
    """
    # Get instances of a data release and SFD dust map.
    decals = DECALS()
    dr = decals.datarelease
    sfd= decals.sfdmap
    cat = dr.catalogue

    mystart = cat.size * comm.rank // comm.size
    myend = cat.size * (comm.rank + 1) // comm.size

    mine = slice(mystart, myend)
    # Define the fluxes, corrected for MW transmission.

    flux  = cat['DECAM_FLUX'][mine].T
    flux /= cat['DECAM_MW_TRANSMISSION'][mine].T
    # Now do the selection ...
    pmask = cat['BRICK_PRIMARY'][mine] == 1

    mask  = cuts.Fluxes.ELG(gflux=flux[1],rflux=flux[2],zflux=flux[4])

    mask &= pmask[None, :]

    select_fraction = sum(comm.allgather(mask.sum(axis=1))) / sum(comm.allgather(pmask.sum()))

    if comm.rank == 0:
        print ('Selected Fraction by Fluxes cuts')

        print ('\n'.join([
            '%s : %g' % v for v in
            zip(cuts.Fluxes.ELG, select_fraction)]))

    mask = mask.all(axis=0)

    total_elg = sum(comm.allgather(mask.sum()))
    total_primary = sum(comm.allgather(pmask.sum()))
    select_fraction = 1.0 * total_elg / total_primary
    if comm.rank == 0:
        print ('Total %d out of %d, ratio=%g' % (total_elg, total_primary, select_fraction))

    # ... and extract only the objects which passed the cuts.
    # At this point we convert fluxes to (extinction corrected)
    # magnitudes, ignoring errors.
    RA   = cat[ 'RA'][mine][mask]
    DEC   = cat['DEC'][mine][mask]
    MAG  = 22.5-2.5*N.log10((flux[:,mask]).clip(1e-15,1e15) )

    # Now we need to pass this through our mask since galaxies can
    # appear even in regions where our nominal depth is insufficient
    # for a complete sample.

    lim = cuts.findlim(dr, sfd, 
                (RA, DEC), 
                ['g','r','z'])

    # the ordering is the same as the call to findlim
    glim,rlim,zlim = lim

    missing_depth = sum(comm.allgather(N.isinf(glim).sum(axis=-1)))

    if comm.rank == 0:
        print('Objects in missing depth images: ', missing_depth)

    mask = cuts.Completeness.ELG(glim=glim,rlim=rlim,zlim=zlim)

    selected_fraction = 1.0 * sum(comm.allgather(mask.sum(axis=1))) \
            / sum(comm.allgather(len(mask.T)))
    if comm.rank == 0:
        print ('Selected Fraction by Completeness cuts')
        print ('\n'.join([
            '%s : %g' % v for v in
            zip(cuts.Completeness.ELG, selected_fraction)]))

    mask = mask.all(axis=0)
    total_complete = sum(comm.allgather(mask.sum()))
    if comm.rank == 0:
        print("total number of complete ELGs:", total_complete)

    RA   = RA [mask]
    DEC   = DEC[mask]
    MAG  = MAG[:, mask]

    RA = N.concatenate(comm.allgather(RA))
    DEC = N.concatenate(comm.allgather(DEC))
    MAG = N.concatenate(comm.allgather(MAG), axis=-1)
    return( (RA,DEC,MAG) )
    #


def diagnostic_plots(ra,dec,mag):
    """
    Makes some "useful" diagnostic plots.
    """
    from matplotlib.figure import Figure
    from matplotlib.backends.backend_agg import FigureCanvasAgg
    g = mag[1, :]
    r = mag[2, :]
    z = mag[4, :]
    fig = Figure()
    ax = fig.add_subplot(111)
    ax.hist2d(ra, dec)
    canvas = FigureCanvasAgg(fig)
    fig.savefig('elg-ra-dec.png')
    #
    fig = Figure()
    ax = fig.add_subplot(111)
    ax.plot(r, r - z, '. ')
    ax.axvline(23.4, label='r < 0.3')
    ax.axhline(0.3, label='(r-z) > 0.3')
    ax.axhline(1.5, label='(r-z) < 1.5')
    ax.set_xlabel('Intrinsic z')
    ax.set_ylabel('Intrinsic r-z')
    ax.grid()
    canvas = FigureCanvasAgg(fig)
    fig.savefig('elg-ri-rz.png')
    #
    fig = Figure()
    ax = fig.add_subplot(111)
    x = N.linspace(
            (g - r).min(),
            (g - r).max(), 100)
    ax.plot(g - r, r - z, '. ')
    ax.plot(x, x + 0.2, label='(g-r) < (r-z) - 0.2')
    ax.plot(x, 1.2 - x , label='(g-r) > 1.2 - (r-z)')
    ax.set_xlabel('Intrinsic g-r')
    ax.set_ylabel('Intrinsic r-z')
    ax.grid()
    canvas = FigureCanvasAgg(fig)
    fig.savefig('elg-gr-rz.png')
    #




if __name__=="__main__":
    from sys import argv
    #
    ra,dc,mag = select_elgs()
    if MPI.COMM_WORLD.rank == 0:
        if len(argv) > 1 and argv[1] == '--plot':
            diagnostic_plots(ra, dc, mag)
        else:
            # Just write the sample to an ascii text file.
            ff = open("elgs.rdz","w")
            ff.write("# %13s %15s %15s %15s %15s %15s\n"%\
              ("RA","DEC","PhotoZ","g","r","z"))
            for i in range(ra.size):
                ff.write("%15.10f %15.10f %15.10f %15.10f %15.10f %15.10f\n"%\
                  (ra[i],dc[i],0.5,mag[1,i],mag[2,i],mag[4,i]))
            ff.close()
        #
