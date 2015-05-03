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
from imaginglss.utils       import sharedmem
from imaginglss             import SFDMap
from imaginglss             import DataRelease
from imaginglss.analysis    import cuts





def select_elgs():
    """
    select_elgs()
    Does the actual selection, imposing cuts on the fluxes
    """
    # Get instances of a data release and SFD dust map.
    dr = DataRelease()
    sfd= SFDMap()
    # Define the fluxes, corrected for MW transmission.
    brickname = dr.catalogue['BRICKNAME']
    flux  = dr.catalogue['DECAM_FLUX'].T
    trn   = dr.catalogue['DECAM_MW_TRANSMISSION'].T
    GFLUX = flux[1] / trn[1]
    RFLUX = flux[2] / trn[2]
    ZFLUX = flux[4] / trn[4]
    # Now do the selection ...
    primary= dr.catalogue['BRICK_PRIMARY']
    pmask = primary == 1

    mask  = cuts.Fluxes.ELG(gflux=GFLUX,rflux=RFLUX,zflux=ZFLUX)
    mask &= pmask[None, :]

    print ('Selected Fraction by Fluxes cuts')

    print ('\n'.join([
        '%s : %g' % v for v in
        zip(cuts.Fluxes.ELG, 1.0 * mask.sum(axis=1) / pmask.sum())]))

    mask = mask.all(axis=0)

    print ('Total %d out of %d, ratio=%g' % (mask.sum(), pmask.sum(), 1.0 * mask.sum() / pmask.sum()))

    # ... and extract only the objects which passed the cuts.
    # At this point we convert fluxes to (extinction corrected)
    # magnitudes, ignoring errors.
    ra   = dr.catalogue[ 'RA'][mask]
    dc   = dr.catalogue['DEC'][mask]
    print ('working on ', len(ra), 'items')
    mag  = 22.5-2.5*N.log10( (flux[:,mask]/trn[:,mask]).clip(1e-15,1e15) )
    # Now we need to pass this through our mask since galaxies can
    # appear even in regions where our nominal depth is insufficient
    # for a complete sample.
    print(N.unique(brickname[mask][:100]))
#    print(dr.images['depth']['z'].open(dr.brickindex.get_brick(325914)))
    with sharedmem.MapReduce() as pool:
        (RA,DEC),arg = dr.brickindex.optimize((ra,dc),return_index=True)
        chunksize = 1024
        def work(i):
            print('', end='.')
#            print(i, '/', len(RA))
            coord = (RA[i:i+chunksize],DEC[i:i+chunksize])
            lim = cuts.findlim(dr,sfd,coord,['g','r','z'])
#            print('done',i)
            return(lim)
        # the ordering is the same as the call to findlim
        glim,rlim,zlim = N.concatenate(\
            pool.map(work,range(0,len(RA),chunksize)),axis=-1)

        print('', end='\n')
        print('glim', glim)
        print('rlim', rlim)
        print('zlim', zlim)

        mask = cuts.Completeness.ELG(glim=glim,rlim=rlim,zlim=zlim)
        print ('Selected Fraction by Completeness cuts')

        print ('\n'.join([
            '%s : %g' % v for v in
            zip(cuts.Completeness.ELG, 1.0 * mask.sum(axis=1) / len(mask))]))

        mask = mask.all(axis=0)
        print(mask.sum())
        ra   = RA [mask]
        dc   = DEC[mask]
        mag  = mag[:,arg[mask]]
    return( (ra,dc,mag) )
    #




def diagnostic_plots(ra,dec,mag):
    """
    diagnostic_plots(ra,dec,mag):
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
