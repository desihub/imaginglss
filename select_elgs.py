#!/usr/bin/env python
#
# Code to select ELG targets from Tractor catalogs.
# The target selection criteria are described at:
#   https://desi.lbl.gov/trac/wiki/TargetSelection
#
from __future__ import print_function


__author__ = "Yu Feng and Martin White"
__version__ = "1.0"
__email__  = "yfeng1@berkeley.edu or mwhite@berkeley.edu"



import numpy as N
from model.datarelease import DataRelease



def select_elgs():
    """
    select_elgs():
    Does the actual selection, imposing cuts on the fluxes
    """
    # Get the catalogs.
    dr = DataRelease()
    # Define the fluxes.
    flux  = dr.catalogue['DECAM_FLUX'].T
    ext   = dr.catalogue['EXTINCTION'].T
    GFLUX = flux[1] * 10 ** (ext[1] / 2.5)
    RFLUX = flux[2] * 10 ** (ext[2] / 2.5)
    ZFLUX = flux[4] * 10 ** (ext[4] / 2.5)
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
    mag  = 22.5-2.5*N.log10(flux[:,mask].clip(1e-15,1e15))-ext[:,mask]
    return( (ra,dc,mag) )
    #


if __name__=="__main__":
    ra,dc,mag = select_elgs()
    print("# %13s %15s %15s %15s %15s"%("RA","DEC","g","r","z"))
    for i in range(ra.size):
        print("%15.10f %15.10f %15.10f %15.10f %15.10f"%\
          (ra[i],dc[i],mag[1,i],mag[2,i],mag[4,i]))
    #
