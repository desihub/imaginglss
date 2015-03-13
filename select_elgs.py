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
# usage: python select_elg.py [--plot]
#
from __future__ import print_function


__author__ = "Yu Feng and Martin White"
__version__ = "1.0"
__email__  = "yfeng1@berkeley.edu or mjwhite@lbl.gov"



import numpy as N
from model.sfdmap      import SFDMap
from model.datarelease import DataRelease



def select_elgs():
    """
    select_elgs():
    Does the actual selection, imposing cuts on the fluxes
    """
    # Get instances of a data release and SFD dust map.
    dr = DataRelease()
    sfd= SFDMap(dustdir="/project/projectdirs/desi/software/edison/dust/v0_0/")
    # Define the fluxes.
    flux  = dr.catalogue['DECAM_FLUX'].T
    trn   = dr.catalogue['DECAM_MW_TRANSMISSION'].T
    GFLUX = flux[1] / trn[1]
    RFLUX = flux[2] / trn[2]
    ZFLUX = flux[4] / trn[4]
    # Now do the selection
    primary = dr.catalogue['BRICK_PRIMARY']
    mask  = (primary == 1)
    mask &= RFLUX > 10**((22.5-23.4) / 2.5)
    mask &= ZFLUX > 10**((0.3) / 2.5) * RFLUX
    mask &= ZFLUX < 10**((1.5) / 2.5) * RFLUX
    mask &= RFLUX ** 2< GFLUX * ZFLUX * 10 ** (-0.2/2.5)
    mask &= ZFLUX > GFLUX * 10**(1.2/2.5)
    # and extract only the objects which passed the cuts.
    # At this point we convert fluxes to (extinction corrected)
    # magnitudes, ignoring errors.
    ra   = dr.catalogue[ 'RA'][mask]
    dc   = dr.catalogue['DEC'][mask]
    mag  = 22.5-2.5*N.log10( (flux[:,mask]/trn[:,mask]).clip(1e-15,1e15) )
    # Now we need to pass this through our mask since galaxies can
    # appear even in regions where our nominal depth is insufficient
    # for a complete sample.  Like in make_random this should probably
    # call out to another routine.  For now I'll just indent it helpfully.
    if True:
        (RA,DEC),arg = dr.brickindex.optimize((ra,dc),return_index=True)
        coord = (RA,DEC)
        # Get the flux depths and MW transmission in the relevant filters.
        # For out-of-bounds points, return 0.
        ebv  = sfd.extinction(None,RA,DEC,get_ebv=True)
        rdep = dr.readout(coord,dr.images['depth']['r'],default=0)
        rtrn = 10.0**(-ebv*dr.extinction['r']/2.5)
        # For now we use a 5-sigma cut in extinction-correct flux as our limit.
        # Recall "depth" is stored as inverse variance.
        rlim = 5.0/N.sqrt(rdep+1e-30) / rtrn
        ww   = N.nonzero( rlim<10.0**((22.5-23.40)/2.5) )[0]
        ra   = RA[ ww]
        dc   = DEC[ww]
        mag  = mag[:,arg[ww]]
    return( (ra,dc,mag) )
    #




def diagnostic_plots(ra, dec, mag):
    from matplotlib.figure import Figure
    from matplotlib.backends.backend_agg import FigureCanvasAgg

    g = mag[:, 1]
    r = mag[:, 2]
    z = mag[:, 4]
    fig = Figure()
    ax = fig.add_subplot(111)
    ax.hist2d(ra, dec)
    canvas = FigureCanvasAgg(fig)
    fig.savefig('elg-ra-dec.png')

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




if __name__=="__main__":
    from sys import argv
    #
    ra,dc,mag = select_elgs()
    if len(argv) > 1 and argv[1] == '--plot':
        diagnostic_plots(ra, dc, mag)
    else:
        ff = open("elgs.rdz","w")
        ff.write("# %13s %15s %15s %15s %15s %15s\n"%\
          ("RA","DEC","PhotoZ","g","r","z"))
        for i in range(ra.size):
            ff.write("%15.10f %15.10f %15.10f %15.10f %15.10f %15.10f\n"%\
              (ra[i],dc[i],0.5,mag[1,i],mag[2,i],mag[4,i]))
        ff.close()
    #
