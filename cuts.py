# "Helper" routines to centralize the flux and color cuts
# in one place where they can be called from other routines.
# Contains two useful classes:
#	Fluxes --	contains methods which return True when
#			an object's fluxes pass a certain cut.
#	Completeness --	contains methods which return True when	a given
#			depth puts it in the 100% complete sample.
# as well as a routine to compute N-sigma flux limits (findlim).

import numpy as N

__author__ = "Yu Feng and Martin White"
__version__ = "1.0"
__email__  = "yfeng1@berkeley.edu or mjwhite@lbl.gov"




def findlim(dr,sfd,coord,bands,sigma=5.0):
    """
    findlim(dr,sfd,coord,bands,sigma=5.0):
    Get the flux depths (corrected for MW transmission) in the relevant
    filters.  For out-of-bounds points, return 0.
    Currently this implements a simple N-sigma cut.
    """
    assert isinstance(bands,(list,tuple))
    ebv  = sfd.ebv(coord[0],coord[1])
    ret = []
    for band in bands:
        rdep = dr.readout(coord,dr.images['depth'][band],\
                          default=0,ignore_missing=True)
        rtrn = 10.0**(-ebv*dr.extinction[band]/2.5)
        # For now we use an N-sigma cut in extinction-correct flux as our limit.
        # Recall "depth" is stored as inverse variance.
        rlim = sigma / N.sqrt(rdep+1e-30) / rtrn
        ret.append(rlim)
    return(ret)
    #


class Fluxes:
    """
    A collection of helpful (static) methods to check whether an object's
    flux passes a given selection criterion (e.g. LRG, ELG or QSO).
    These cuts assume we are passed the extinction-corrected fluxes
    (flux/mw_transmission) and are taken from:
      https://desi.lbl.gov/trac/wiki/TargetSelection
    The requirement on BRICK_PRIMARY is handled elsewhere, these methods just
    look after the flux requirements.
    """
    @staticmethod
    def LRG(rflux,zflux,wflux):
        mask  = rflux > 10**((22.5-23.00)/2.5)
        mask &= zflux > 10**((22.5-20.56)/2.5)
        mask &= wflux > 10**((22.5-19.50)/2.5)
        mask &= zflux > 10**((1.6)       /2.5)*rflux
        mask &= wflux*rflux**(1.33-1)>zflux**1.33*10**(-0.33/2.5)
        return(mask)
    #
    @staticmethod
    def ELG(gflux,rflux,zflux):
        mask  = rflux    > 10**((22.5-23.4)/2.5)
        mask &= zflux    > 10**((0.3)      /2.5)*rflux
        mask &= zflux    < 10**((1.5)      /2.5)*rflux
        mask &= rflux**2 < gflux*zflux*10**(-0.2/2.5)
        mask &= zflux    > gflux*10**(1.2/2.5)
        return(mask)
    #
    @staticmethod
    def QSO(gflux,rflux,w1flux,w2flux):
        wflux = 0.67*w1flux + 0.33*w2flux
        mask  = rflux > 10**((22.5-23.0)/2.5)
        mask &= rflux < 10**((1.0)      /2.5)*gflux
        mask &= wflux*gflux**1.2 > 10**(2./2.5)*rflux**(1+1.2)
        return(mask)
    @staticmethod
    def BGS(rflux):
        mask  = rflux > 10**((22.5-19.5)/2.5)
        return(mask)
    #
 


class Completeness:
    """
    A collection of helpful (static) methods to check whether a given
    position has imaging deep enough to be 100% complete to a certain
    selection (e.g. LRG, ELG or QSO).
    We should probably modify this later to easily allow querying of
    brighter limits.
    The methods assume they are passed N-sigma limits, corrected for
    MW transmission.
    Note: for a region to pass the N-sigma "noise" must be LESS THAN the
    required depth.
    The requirements on colors are tricky to implement here.
    """
    @staticmethod
    def LRG(rlim,zlim,wlim):
        # This is not finished yet.
        raise UnimplementedError
        mask  = rlim<10**((22.5-23.00    )/2.5)
        mask &= zlim<10**((22.5-20.56    )/2.5)
        mask &= wlim<10**((22.5-19.50    )/2.5)
        mask &= zlim<10**((22.5-23.00+1.6)/2.5)
        return(mask)
    #
    @staticmethod
    def ELG(glim,rlim,zlim):
        mask  = rlim < 10**((22.5-23.4        )/2.5) 
        mask &= zlim < 10**((22.5-23.4+0.3    )/2.5)
        mask &= glim < 10**((22.5-23.4-1.5+0.2)/2.5)
        return(mask)
    #
    @staticmethod
    def QSO(glim,rlim,w1lim,w2lim):
        # This is not finished yet.
        raise UnimplementedError
        mask  = rlim<10**((22.5-23.00    )/2.5)
        mask &= glim<10**((22.5-23.00-1.0)/2.5)
        return(mask)
    #
    @staticmethod
    def BGS(rlim):
        mask  = rlim<10**((22.5-19.5)/2.5)
        return(mask)
    #
