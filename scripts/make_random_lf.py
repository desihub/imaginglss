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

import os.path; import sys; sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..'))

from argparse import ArgumentParser

ap = ArgumentParser("make_random.py")
ap.add_argument("Nran", type=int, help="Minimum number of randoms")
ap.add_argument("ObjectType", choices=["QSO", "LRG", "ELG", "BGS"])
ap.add_argument("output")
ap.add_argument("--conf", default=None, 
        help="Path to the imaginglss config file, default is from DECALS_PY_CONFIG")

ns = ap.parse_args()

import numpy             as np
from   imaginglss             import DECALS
from   imaginglss.analysis    import cuts
from   imaginglss.model.datarelease import Footprint
from   imaginglss.utils import filehandler
from mpi4py            import MPI

np.seterr(divide='ignore', invalid='ignore')

def select_objs(sampl, configfile, comm=MPI.COMM_WORLD):
    """
    Does the actual selection, imposing cuts on the fluxes
    """
    # Get instances of a data release and SFD dust map.
    decals = DECALS(configfile)
    dr     = decals.datarelease
    sfd    = decals.sfdmap
    cat    = dr.catalogue
    #
    mystart = cat.size * comm.rank // comm.size
    myend = cat.size * (comm.rank + 1) // comm.size
    #
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
            '%-50s : %.4f' % v for v in
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
            '%-50s : %.4f' % v for v in
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
        MAG[band] = 22.5-2.5*np.log10((flux[band]).clip(1e-15,1e15) )
    for band in flux:
        flux[band] = np.concatenate(comm.allgather(flux[band]))
        MAG[band] = np.concatenate(comm.allgather(MAG[band]))
    RA  = np.concatenate(comm.allgather(RA ))
    DEC = np.concatenate(comm.allgather(DEC))
    return( (RA,DEC,flux) )

def fill_random(footprint, Nran, rng):
    """
    Generate uniformly distributed points within the boundary that lie in
    bricks.  We generate in the ra/dec area, then remove points not in any
    bricks.  This hugely increases the memory efficiency for footprints,
    like DR1, where the survey covers several disjoint patches of the sky.

    """


    coord = np.empty((2, Nran))

    ramin,ramax,dcmin,dcmax = footprint.range

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


def apply_samp_cut(coord, dr, sfd, samp):
    """
    apply_samp_cut(coord, dr, sfd, samp): 
    """

    return mask, rmag, cut

def sample_lf(lims, bands, pool, rng, Nran=10000):
    """
        Sample randoms from luminosity function.

        Parameters
        ----------
        lims : dict of array
            each array is the limit read out from findlim for that band
            the length of the array is the total number of positions on the sky
            that have been used.
        
        bands : list of string
            bands to compute the completeness 

        pool : dict of array
            each item is the flux of objects to selected. This is a realizatoin
            of the LF of the population to be mocked.

        rng  : :py:class:`numpy.random.RandomState`
            Random number generator

        Nran : int
            Number of random objects to select per position.

        Returns
        -------
        completeness: 
            list (per band) of array of completeness per random point
    """

    indices = rng.randint(len(pool[bands[0]]), size=len(lims[bands[0]]) * Nran)
    completeness = []
    for band in bands:
        expanded_lims = np.repeat(lims[band], Nran)
        flux = pool[band][indices]
        print(band)
        print(expanded_lims)
        print(flux)
        completeness.append(1.0 * (expanded_lims < flux).reshape(-1, Nran).sum(axis=-1) / Nran)
    return np.array(completeness)
        

def make_random(samp, Nran, pool, configfile, comm=MPI.COMM_WORLD):
    """
    Does the work of making randoms.  The sample type is passed as a string.
    """

    # Get the total footprint bounds, to throw randoms within, and an E(B-V)
    # map instance.
    decals = DECALS(configfile)
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

    # Distribute myNran, proportional to the number of bricks
    # important if a rank has no bricks, it expects no randoms

    myNran = rng.poisson(footprint.area / dr.footprint.area * Nran)
    
    print (comm.rank, 'has', len(mybricks), 'bricks', myNran, 'randoms')
    # fill it with random points 
    coord = fill_random(footprint, myNran, rng)

    Nran = sum(comm.allgather(len(coord[0])))

    cut = getattr(cuts.Completeness, samp)

    lim = cuts.findlim(dr,sfd,coord, cut.bands)

    completeness = sample_lf(lim, cut.bands, pool, rng, Nran=100)

    # It's also useful to have r magnitude later.
    rmag = 22.5-2.5*np.log10( lim['r'].clip(1e-15,1e15) )

    completeness = comm.gather(completeness)
    coord = comm.gather(coord)
    rmag = comm.gather(rmag)

    if comm.rank == 0:
        completeness = np.concatenate(completeness, axis=-1)
        completeness = completeness.prod(axis=0)
        coord = np.concatenate(coord, axis=-1)
        rmag = np.concatenate(rmag)
        fraction = len(coord[0]) * 1.0 / Nran
        print('Accept rate', fraction)
        print('Total area (sq.deg.) ',dr.footprint.area * fraction)
        print('Done!')

    return coord, completeness, rmag
    #


if __name__ == '__main__':    

    RA, DEC, fluxes = select_objs(ns.ObjectType, configfile=ns.conf)

    coord, completeness, rmag = make_random(ns.ObjectType, pool=fluxes, configfile=ns.conf, Nran=ns.Nran)

    if MPI.COMM_WORLD.rank == 0:
#        filehandler.write(ns.output, dict(
#                ra=coord[0], 
#                dec=coord[1], 
#                weights=, 

        with open(ns.output,'w') as ff:
            ff.write("# ra dec weight rmag\n")
            for j in range(0, len(rmag)):
                ff.write("%15.10f %15.10f %15.10f %15.10f\n"%\
                  (coord[0][j],coord[1][j], completeness[j],rmag[j]))

