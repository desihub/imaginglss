#!/usr/bin/env python
#
# Code to select targets from Tractor catalogs.
# The target selection criteria are described at:
#   https://desi.lbl.gov/trac/wiki/TargetSelection
#
# These catalogs also need to be vetod using the
# depths and masks which are part of Tractor, since
# catalog entries can appear in places where the
# nominal depth is insufficient to be complete.
#
# Usage: python select_objs.py [--plot]
#
from __future__ import print_function


__author__ = "Yu Feng and Martin White"
__version__ = "1.0"
__email__  = "yfeng1@berkeley.edu or mjwhite@lbl.gov"

import os.path; import sys; sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..'))

import numpy as np
from imaginglss             import DECALS
from imaginglss.analysis    import cuts
from imaginglss.analysis    import targetselection

from mpi4py import MPI

from argparse import ArgumentParser

ap = ArgumentParser("select_objs.py")
ap.add_argument("ObjectType", choices=["QSO", "LRG", "ELG", "BGS"])
ap.add_argument("output")
ap.add_argument("--sigma-z", type=float, default=3.0)
ap.add_argument("--sigma-g", type=float, default=5.0)
ap.add_argument("--sigma-r", type=float, default=5.0)
ap.add_argument("--conf", default=None, 
        help="Path to the imaginglss config file, default is from DECALS_PY_CONFIG")

ns = ap.parse_args()


np.seterr(divide='ignore', invalid='ignore')

def select_objs(sampl, conffile, comm=MPI.COMM_WORLD):
    """
    Does the actual selection, imposing cuts on the fluxes
    """
    # Get instances of a data release and SFD dust map.
    decals = DECALS(conffile)
    dr     = decals.datarelease
    sfd    = decals.sfdmap
    cat    = dr.catalogue
    #
    mystart = cat.size * comm.rank // comm.size
    myend = cat.size * (comm.rank + 1) // comm.size
    #
    mine = slice(mystart, myend)

    with dr.catalogue as cat:
        rows = cat[mine]
        fluxcut = getattr(targetselection, sampl)
        for expr in fluxcut:
            mask = expr.visit(rows)
            selected = sum(comm.allgather(mask.sum())) 
            if comm.rank == 0:
                print("%s : %d / %d = %g"
                    % (
                    str(expr), selected, cat.size, 
                    1. * selected / cat.size))
        mask = fluxcut.visit(rows)
        selected = sum(comm.allgather(mask.sum())) 
        if comm.rank == 0:
            print("%s : %d / %d = %g"
                % (
                str(fluxcut), selected, cat.size, 
                1. * selected / cat.size))

    total_elg = sum(comm.allgather(mask.sum()))
    # ... and extract only the objects which passed the cuts.
    RA   = cat[ 'RA'][mine][mask]
    DEC   = cat['DEC'][mine][mask]
    FLUX = cat['DECAM_FLUX'][mine][mask] / cat['DECAM_MW_TRANSMISSION'][mine][mask]
    MAG = 22.5-2.5*np.log10((FLUX).clip(1e-15,1e15))

    # Now we need to pass this through our mask since galaxies can
    # appear even in regions where our nominal depth is insufficient
    # for a complete sample.
    compcut = getattr(cuts.Completeness, sampl)
    d = dict(z=ns.sigma_z, g=ns.sigma_g, r=ns.sigma_r)
    sigma = [ d[band] for band in compcut.bands]
    lim = cuts.findlim(dr, sfd, 
                (RA, DEC), 
                compcut.bands, sigma=sigma)
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
            '%-50s : %.4f' % v for v in
            zip(compcut, selected_fraction)]))
    mask = mask.all(axis=0)
    total_complete = sum(comm.allgather(mask.sum()))
    if comm.rank == 0:
        print('Total number of objects in complete area: ', total_complete)
    RA  = RA [mask]
    DEC = DEC[mask]
    MAG = MAG[mask]
    RA  = np.concatenate(comm.allgather(RA ))
    DEC = np.concatenate(comm.allgather(DEC))
    MAG = np.concatenate(comm.allgather(MAG), axis=0)
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

    ra,dc,mag = select_objs(ns.ObjectType, ns.conf)

    if MPI.COMM_WORLD.rank == 0:
        # Just write the sample to an ascii text file.
        ff = file(ns.output,"w")
        with ff:
            ff.write("# sigma_z=%g sigma_g=%g sigma_r=%g\n" % (ns.sigma_z, ns.sigma_g, ns.sigma_r))
            ff.write("# %13s %15s %15s" % ("RA","DEC","PhotoZ"))
            for band in 'ugrizY':
                ff.write(" %15s" % band)
            ff.write("\n")
            for i in range(ra.size):
                ff.write("%15.10f %15.10f %15.10f"
                    % (ra[i],dc[i],0.5))
                for band in range(6):
                    ff.write(" %15.10f" % mag[i][band])
                ff.write("\n")
        #
