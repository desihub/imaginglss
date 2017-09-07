"""
    These cuts assume we are passed the extinction-corrected fluxes
    (flux/mw_transmission) and are taken from:

      https://desi.lbl.gov/trac/wiki/TargetSelection
"""

__author__ = "Yu Feng and Martin White"
__version__ = "1.0"
__email__  = "yfeng1@berkeley.edu or mjwhite@lbl.gov"

WFLUX = 0.75 * W1FLUX + 0.25 * W2FLUX
GRZFLUX = (GFLUX + 0.8* RFLUX + 0.5* ZFLUX ) / 2.4
SNRW1 = (W1FLUX * W1FLUX_IVAR ** 0.5)
SNRW2 = (W2FLUX * W2FLUX_IVAR ** 0.5)

LRG =  BRICK_PRIMARY != 0
LRG &= ZFLUX > 10**((22.5-20.46)/2.5)
LRG &= ZFLUX > RFLUX * 10**(1.5/2.5)
LRG &= W1FLUX * RFLUX ** (1.8-1) > ZFLUX**1.8 * 10**(-1.0/2.5)
LRG &= W1FLUX > 0

ELG =  BRICK_PRIMARY != 0
ELG &= RFLUX > 10**((22.5-23.4)/2.5)
ELG &= ZFLUX > 10**(0.3/2.5) * RFLUX
ELG &= ZFLUX < 10**(1.6/2.5) * RFLUX
ELG &= RFLUX**2.15 < GFLUX * ZFLUX**1.15 * 10**(-0.15/2.5)
ELG &= ZFLUX**1.2 < GFLUX * RFLUX**0.2 * 10**(1.6/2.5)
#ELG &= Max(SHAPEDEV_R, SHAPEEXP_R) < 1.5

# QSO by colors only
QSOC  = RFLUX > 10**((22.5-23.0)/2.5)
QSOC &= GRZFLUX < 10**((22.5-17.0)/2.5)
QSOC &= RFLUX < 10**(1.3/2.5) * GFLUX
QSOC &= ZFLUX > 10**(-0.3/2.5) * RFLUX
QSOC &= ZFLUX < 10**(1.1/2.5) * RFLUX
QSOC &= WFLUX * GFLUX > 10**(-1.0/2.5) * ZFLUX * GRZFLUX
QSOC &= W2FLUX > W1FLUX * 10**(-0.4 / 2.5)
QSOC &= SNRW1 > 4
QSOC &= SNRW2 > 2

QSO = BRICK_PRIMARY != 0
QSO &= QSOC
QSO &= TYPE == 'PSF '

# David's variant of QSO
QSOd = BRICK_PRIMARY != 0
QSOd &= QSOC
QSOd &= Max(SHAPEDEV_R, SHAPEEXP_R) < 0.5

BGS =  BRICK_PRIMARY != 0
BGS &= TYPE != 'PSF '
BGS &= RFLUX > 10**((22.5-19.5)/2.5)


