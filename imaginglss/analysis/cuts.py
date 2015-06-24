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
    Get the flux depths (corrected for MW transmission) in the relevant
    filters.  

    For out-of-bounds points, return 0.
    Currently this implements a simple N-sigma cut.

    Parameters
    ----------
    dr : :py:class:`model.datarelease.DataRelease`
        data release
    sfd  : :py:class:`model.sfdmap.SFDMap`
        dust extinction map
    coord :  array_like
        coord = (RA, DEC)
    bands : list
        ['g', 'u', 'z', ...]
    sigma : float, or list
        confidence level for each band

    Returns
    -------
    lims : dict 
        one lim per band.

    """
    assert isinstance(bands,(list,tuple))
    ebv  = sfd.ebv(coord[0],coord[1])
    ret = {}
    if N.isscalar(sigma):
        sigma = [sigma] * len(bands)
    for s, band in zip(sigma, bands):
        rtrn = 10.0**(-ebv*dr.extinction[band]/2.5)

        rdep = dr.readout(coord,dr.images['depth'][band],\
                          default=+0.0,ignore_missing=True)

        # For now we use an N-sigma cut in extinction-correct flux as our limit.
        # Recall "depth" is stored as inverse variance.
        
        # Bricks with missing depth images are simply passed as 0.0
        # which would become +inf
        # we will use this to detect objects on missing bricks.

        rlim = s * rdep ** -0.5 / rtrn
        ret[band] = rlim
    return(ret)
    #

class Cuts(object):
    """ A convenient representation of a set of Cuts.
        
        For the ease of self-introspection, a :py:class:`Cuts` object 
        stores the cut expressions as a list of strings. 
        Python provides this useful function :code:`eval` that evaluates the
        expression with a set of variables in a local scope. We make use
        of :code:`eval` to evaluate these expressions for each cut criteria. 

        These variables used in the cuts evaluation are decleared in 'bands' parameter.
        Some times a cut criteria needs a band flux that is constructed from other
        band fluxes, these can be provided in 'transforms' parameter.

        Attributes
        ----------
        cuts : list
            a list of python expressions as strings. 
            For example,  ["r > 10**((22.5-23.00)/2.5)"]

        transforms : kwargs
            a dict of band = "expression"            
            where "expressions" is like "w = 0.75*w1 + 0.25*w2",

        bands : list
            a list of bands that this Cuts object uses.

        Parameters
        ----------
        **local : dict
            the arguments to eval the cut expression. for example, rflux=[1.1, 2.3, ...]
        
        Examples
        --------
        
        The following code defines a QSO color cut by r, g, and W band fluxes,
        where the W band flux is obtained by interpolating W1, and W2 band fluxes.

        >>> QSO = Cuts([
        "r > 10**((22.5-23.0)/2.5)",
        "r < 10**((1.0)      /2.5)*g",
        "W*g**1.2 > 10**(-0.4/2.5)*r**(1+1.2)", 
        ], bands=['r', 'g', 'W1', 'W2'],
        W="0.75*W1 + 0.25*W2",
        )

        >>> mask = QSO(r=rflux, g=gflux, W1=W1flux, W2=W2flux)
        >>> mask.shape
        (3, ....)
        
        >>> for m, cut in zip(mask, QSO):
        >>>    print cut, ' : ', m.sum()
        r > 10**((22.5-23.0)/2.5) :  300
        r < 10**((1.0)      /2.5)*g : 400
        W*g**1.2 > 10**(-0.4/2.5)*r**(1+1.2)  : 100
          
        Returns
        -------
        mask : array_like (Ncuts, Nitem)
            True if item passed the cut, False if not. 

    """
    def __init__(self, cuts, bands, **transforms):
        self.cuts = cuts
        self.transforms = transforms
        self.bands = bands
        self.__doc__ = self.__class__.__doc__ + \
        """
        Notes
        -----
        Bands are : %(bands)s
        """ % dict(bands= str(self.bands))

    def __len__(self):
        return len(self.cuts)

    def __call__(self, **fluxes):

        # evaluate the cuts for the input fluxes as keyword arguments
        
        # first build the transformed flux variables (eg. W from W1 and W2 for 
        # QSOs

        for val in self.transforms:
            fluxes[val] = eval(self.transforms[val], fluxes)

        for band in self.bands:
            if band not in fluxes:
                raise ValueError("band '%s' is not specified" % band)
 
        return N.array([
            eval(expr, fluxes)
            for expr in self.cuts])

    def __iter__(self):
        return iter(self.cuts)

    def __repr__(self):
        return "Cut in %s\n%s\n%s" % (repr(self.bands),
                '\n'.join(['%s = %s' % pair 
                    for pair in self.transforms.items()]), 
                '\n'.join(self.cuts))
    
from imaginglss.model.catalogue import C

class ObjectTypes:
    """
    A collection of helpful (static) methods to check whether an object's
    flux passes a given selection criterion (e.g. LRG, ELG or QSO).
    These cuts assume we are passed the extinction-corrected fluxes
    (flux/mw_transmission) and are taken from:
      https://desi.lbl.gov/trac/wiki/TargetSelection
    The requirement on BRICK_PRIMARY is handled elsewhere, these methods just
    look after the flux requirements.
    """
    DECAM_FLUX = C('DECAM_FLUX').T
    DECAM_MW_TRANSMISSION = C('DECAM_MW_TRANSMISSION').T
    WISE_FLUX = C('WISE_FLUX').T
    WISE_MW_TRANSMISSION = C('WISE_MW_TRANSMISSION').T
    BRICK_PRIMARY = C('BRICK_PRIMARY')
    TYPE = C('TYPE')

    GFLUX = DECAM_FLUX[1] / DECAM_MW_TRANSMISSION[1]
    RFLUX = DECAM_FLUX[2] / DECAM_MW_TRANSMISSION[2]
    ZFLUX = DECAM_FLUX[4] / DECAM_MW_TRANSMISSION[4]
    W1FLUX = WISE_FLUX[0] / WISE_MW_TRANSMISSION[0] 
    WFLUX = 0.75 * WISE_FLUX[0] / WISE_MW_TRANSMISSION[0] \
          + 0.25 * WISE_FLUX[1] / WISE_MW_TRANSMISSION[1] 

    LRG =  BRICK_PRIMARY != 0
    LRG &= RFLUX > 10**((22.5-23.0)/2.5)
    LRG &= ZFLUX > 10**((22.5-20.56)/2.5)
    LRG &= W1FLUX > 10**((22.5-19.35)/2.5)
    LRG &= ZFLUX > RFLUX * 10**(1.6/2.5)
    LRG &= W1FLUX * RFLUX ** (1.33-1) > ZFLUX**1.33 * 10**(-0.33/2.5)

    ELG =  BRICK_PRIMARY != 0
    ELG &= RFLUX > 10**((22.5-23.4)/2.5)
    ELG &= ZFLUX > 10**(0.3/2.5) * RFLUX
    ELG &= ZFLUX < 10**(1.5/2.5) * RFLUX
    ELG &= RFLUX**2 < GFLUX * ZFLUX * 10**(-0.2/2.5)
    ELG &= ZFLUX < GFLUX * 10**(1.2/2.5) 

    QSO =  BRICK_PRIMARY != 0
    QSO &= RFLUX > 10**((22.5-23.0)/2.5)
    QSO &= RFLUX < 10**(1.0/2.5) * GFLUX
    QSO &= ZFLUX > 10**(-0.3/2.5) * RFLUX
    QSO &= ZFLUX < 10**(1.1/2.5) * RFLUX
    QSO &= WFLUX * GFLUX**1.2 > 10**(2/2.5) * RFLUX**(1+1.2)

    BGS =  BRICK_PRIMARY != 0
    BGS &= TYPE != 'PSF'
    BGS &=  RFLUX > 10**((22.5-19.35)/2.5)

class Fluxes:
    LRG = Cuts([
        "r > 10**((22.5-23.00)/2.5)",
        "z > 10**((22.5-20.56)/2.5)",
        "W1 > 10**((22.5-19.50)/2.5)",
        "z > 10**((1.6)       /2.5)*r",
        "W1 *r **(1.33-1)>z **1.33*10**(-0.33/2.5)",
        ], bands=['r', 'z', 'W1'])

    SEL1 = Cuts([
        "z > 10**((22.5-19.2)/2.5)",
        "r > z * 10**((-1.7)/2.5)",
        "r < z * 10**((-0.9)/2.5)",
        "z < W1 * 10**((-0.5)/2.5)",
        "z > W1 * 10**((-1.3)/2.5)",
        "z*z < r*W1 * 10**((0.55)/2.5)",
        "z*z > r*W1 * 10**((0.1)/2.5)",
        ], bands=['r', 'z', 'W1'])

    ELG = Cuts([
        "r > 10**((22.5-23.4)/2.5)",
        "z > 10**((0.3)      /2.5)*r",
        "z < 10**((1.5)      /2.5)*r",
        "r **2 < g *z *10**(-0.2/2.5)",
        "z < g *10**(1.2/2.5)",
        ], bands=['r', 'z', 'g'])

    QSO = Cuts([
        "r > 10**((22.5-23.0)/2.5)",
        "r < 10**((1.0)      /2.5)*g",
        "W*g**1.2 > 10**(-0.4/2.5)*r**(1+1.2)", 
        ], bands=['r', 'g', 'W1', 'W2'],
        W="0.75*W1 + 0.25*W2",
        )

    BGS = Cuts([
        "r > 10**((22.5-19.5)/2.5)",
        ], bands = ['r']
        )
 


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
    The requirements on colors are tricky to implement here. We currently
    require a 100% completeness.
    
    """

    LRG = Cuts([
        "r<10**((22.5-23.00    )/2.5)",
        "z<10**((22.5-20.56    )/2.5)",
        "z<10**((22.5-23.00+1.6)/2.5)",
        #"W<10**((22.5-19.50    )/2.5)",
        # "raise UnimplementedError",
        ], bands=['r', 'z'])

    ELG = Cuts([
        "g < 10**((22.5-23.4-1.5+0.2)/2.5)",
        "r < 10**((22.5-23.4        )/2.5)",
        "z < 10**((22.5-23.4+0.3    )/2.5)",
        ], bands=['g', 'z', 'r'])

    QSO = Cuts([
        #"raise UnimplementedError",
        "r<10**((22.5-23.00    )/2.5)",
        "g<10**((22.5-23.00-1.0)/2.5)",
        ], bands=['r', 'g'])
    BGS = Cuts([
        "r<10**((22.5-19.5)/2.5)",
        ], bands=['r'])

    SEL1 = Cuts([
        "r < 10**((22.5-19.2 - 1.7)/2.5)",
        "z < 10**((22.5-19.2)/2.5)",
        ], bands=['r', 'z'])
    #
