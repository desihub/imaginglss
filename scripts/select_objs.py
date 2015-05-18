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

import numpy as np
from imaginglss             import DECALS
from imaginglss.analysis    import cuts

from mpi4py import MPI

from argparse import ArgumentParser

ap = ArgumentParser("select_obj.py")
ap.add_argument("ObjectType", choices=["QSO", "LRG", "ELG", "BGS"])
ap.add_argument("output")
ap.add_argument("--conf", default=None, 
        help="Path to the imaginglss config file, default is from DECALS_PY_CONFIG")

ns = ap.parse_args()


np.seterr(divide='ignore', invalid='ignore')

def select_elgs(sampl, conffile, comm=MPI.COMM_WORLD):
    """
    Does the actual selection, imposing cuts on the fluxes
    """
    # Get instances of a data release and SFD dust map.
    decals = DECALS(conffile)
    dr = decals.datarelease
    sfd= decals.sfdmap
    cat = dr.catalogue

    mystart = cat.size * comm.rank // comm.size
    myend = cat.size * (comm.rank + 1) // comm.size

    mine = slice(mystart, myend)
    # Define the fluxes, corrected for MW transmission.
    decam_flux = (cat['DECAM_FLUX'][mine]/cat['DECAM_MW_TRANSMISSION'][mine]).T
    wise_flux  = (cat[ 'WISE_FLUX'][mine]/cat[ 'WISE_MW_TRANSMISSION'][mine]).T
    # Now do the selection ...
    pmask   = cat['BRICK_PRIMARY'][mine] == 1
    fluxcut = getattr(cuts.Fluxes, sampl)
    flux    = {}
    for band in fluxcut.bands:
        if band in 'ugrizY':
            flux[band] = decam_flux[dr.bands[band]]
        elif band in ['W1', 'W2', 'W3', 'W4']:
            flux[band] = wise_flux[int(band[1]) - 1]
    mask  = fluxcut(**flux)
    mask &= pmask[None, :]
    select_fraction = (1.0*sum(comm.allgather(mask.sum(axis=1)))) / sum(comm.allgather(pmask.sum()))
    if comm.rank == 0:
        print('Selected fraction by flux cuts:')
        print('\n'.join([
            '%-50s : %.3f' % v for v in
            zip(fluxcut, select_fraction)]))
    mask = mask.all(axis=0)
    total_elg = sum(comm.allgather(mask.sum()))
    total_primary = sum(comm.allgather(pmask.sum()))
    select_fraction = 1.0 * total_elg / total_primary
    if comm.rank == 0:
        print('Selected %d out of %d objects (%g)' %\
              (total_elg, total_primary, select_fraction))
    # ... and extract only the objects which passed the cuts.
    # At this point we convert fluxes to (extinction corrected)
    # magnitudes, ignoring errors.
    RA   = cat[ 'RA'][mine][mask]
    DEC   = cat['DEC'][mine][mask]
    for band in flux:
        flux[band] = flux[band][mask]
    # Now we need to pass this through our mask since galaxies can
    # appear even in regions where our nominal depth is insufficient
    # for a complete sample.
    compcut = getattr(cuts.Completeness, sampl)
    lim = cuts.findlim(dr, sfd, 
                (RA, DEC), 
                compcut.bands)
    for band in lim:
        missing_depth = sum(comm.allgather(np.isinf(lim[band]).sum()))
        if comm.rank == 0:
            print('Objects in bricks with missing depth images (',band,'): ',\
                  missing_depth)
    mask = compcut(**lim)
    selected_fraction = 1.0 * sum(comm.allgather(mask.sum(axis=1)))\
            / sum(comm.allgather(len(mask.T)))
    if comm.rank == 0:
        print ('Selected fraction by completeness cuts:')
        print ('\n'.join([
            '%s : %.3f' % v for v in
            zip(compcut, selected_fraction)]))
    mask = mask.all(axis=0)
    total_complete = sum(comm.allgather(mask.sum()))
    if comm.rank == 0:
        print('Total number of objects in complete area: ', total_complete)
    RA  = RA [mask]
    DEC = DEC[mask]
    MAG = {}
    for band in flux:
        flux[band] = flux[band][mask]
    for band in flux:
        MAG[band] = 22.5-2.5*np.log10((flux[band]).clip(1e-15,1e15) )
        MAG[band] = np.concatenate(comm.allgather(MAG[band]))
    RA  = np.concatenate(comm.allgather(RA ))
    DEC = np.concatenate(comm.allgather(DEC))
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
    x = np.linspace(
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

    ra,dc,mag = select_elgs(ns.ObjectType, ns.conf)

    if MPI.COMM_WORLD.rank == 0:
        # Just write the sample to an ascii text file.
        ff = file(ns.output,"w")
        with ff:
            ff.write("# %13s %15s %15s" % ("RA","DEC","PhotoZ"))
            for band in mag:
                ff.write(" %15s" % band)
            ff.write("\n")
            for i in range(ra.size):
                ff.write("%15.10f %15.10f %15.10f"
                    % (ra[i],dc[i],0.5))
                for band in mag:
                    ff.write(" %15.10f" % mag[band][i])
                ff.write("\n")
        #
