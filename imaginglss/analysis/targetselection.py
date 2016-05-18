"""
    Target Selection

    https://desi.lbl.gov/trac/wiki/TargetSelection

    A collection of helpful (static) methods to check whether an object's
    flux passes a given selection criterion (e.g. LRG, ELG or QSO).

    These cuts assume we are passed the extinction-corrected fluxes
    (flux/mw_transmission) and are taken from:

      https://desi.lbl.gov/trac/wiki/TargetSelection

    Put your own target selection definitions in

        local-targetselection.py 

    of the same path.

    Remember to append the name of the object type to __all__ variable
    with __all__.append("ObjectType")
"""

__author__ = "Yu Feng and Martin White"
__version__ = "1.0"
__email__  = "yfeng1@berkeley.edu or mjwhite@lbl.gov"

from imaginglss.model.columnnames import *
from imaginglss.utils.npyquery import Max, Min

LRG =  BRICK_PRIMARY != 0
LRG &= ZFLUX > 10**((22.5-20.46)/2.5)
LRG &= ZFLUX > RFLUX * 10**(1.5/2.5)
LRG &= W1FLUX * RFLUX ** (1.8-1) > ZFLUX**1.8 * 10**(-1.0/2.5)
LRG &= W1FLUX > 0
LRG.limit_bands = 'z'
LRG.bands = 'zr'

ELG =  BRICK_PRIMARY != 0
ELG &= RFLUX > 10**((22.5-23.4)/2.5)
ELG &= ZFLUX > 10**(0.3/2.5) * RFLUX
ELG &= ZFLUX < 10**(1.6/2.5) * RFLUX
ELG &= RFLUX**2.15 < GFLUX * ZFLUX**1.15 * 10**(-0.15/2.5)
ELG &= ZFLUX**1.2 < GFLUX * RFLUX**0.2 * 10**(1.6/2.5)
#ELG &= Max(SHAPEDEV_R, SHAPEEXP_R) < 1.5
ELG.limit_bands = 'r'
ELG.bands = 'rgz'

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
QSO.limit_bands = 'r'
QSO.bands = 'rgz'

# David's variant of QSO
QSOd = BRICK_PRIMARY != 0
QSOd &= QSOC
QSOd &= Max(SHAPEDEV_R, SHAPEEXP_R) < 0.5
QSOd.limit_bands = 'r'
QSOd.bands = 'rgz'

BGS =  BRICK_PRIMARY != 0
BGS &= TYPE != 'PSF '
BGS &= RFLUX > 10**((22.5-19.5)/2.5)
BGS.limit_bands = 'r'
BGS.bands = 'r'


__all__ = []

def load(path):
    """ Recursively load target definitions in the path."""
    import os
    if os.path.isdir(path):
        for f in os.listdir(path):
            if os.path.isdir(path) or \
              (f.endswith('.py') and not f.startswith('_')):
                load(os.path.join(path, f))
    else:
        script = open(path, 'r').read()

        exec(compile(script, path, 'exec'), globals())
        _prune()

# now we try to import a local version the file
def _local():

    import os

    local = os.path.join(os.path.dirname(__file__), 'local-%s' % os.path.basename(__file__))
    if local.endswith('.pyc'): local = local[:-1]

    if os.path.exists(local):
        load(local)

def _prune():
    global __all__
    g = globals()

    # This will filter out names that does not appear to be a target type.
    import imaginglss.model.columnnames as columnnames
    blacklist = dir(columnnames)
    blacklist.append('load')

    __all__ = []

    for k in g.keys():
        if k in blacklist:
            continue
        if k in ['Min', 'Max']:
            continue
        if k.startswith('_'):
            continue
        __all__.append(k)


_local()
_prune()
