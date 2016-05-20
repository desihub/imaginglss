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

# * is evil but exactly what we want to do here
from imaginglss.model.columnnames import *

WFLUX = 0.75 * W1FLUX + 0.25 * W2FLUX[1]
GRZFLUX = (GFLUX + 0.8* RFLUX + 0.5* ZFLUX ) / 2.4
SNRW1 = (WISE_FLUX[0] * WISE_FLUX_IVAR[0] ** 0.5)
SNRW2 = (WISE_FLUX[1] * WISE_FLUX_IVAR[1] ** 0.5)

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


__all__ = []

def load(path):
    """ Recursively load target definitions in the path."""
    import os
    if path.startswith('_'): return
    if os.path.isdir(path):
        for f in sorted(os.listdir(path)):
            load(os.path.join(path, f))
    else:
        if not path.endswith('.py'): return
        script = open(path, 'r').read()
        d = {}
        exec(compile(script, path, 'exec'), globals(), d)
        _export(d)

def _gather_color_bands(expr):
    result = []
    def walk(expr):
        for c in expr.children:
            walk(c)
        if hasattr(expr, 'band'):
            result.append(expr.band)
    walk(expr)
    return list(set(result))

def _gather_magnitude_bands(expr):
    from imaginglss.utils.npyquery import Expr, Literal
    result = []
    def walk(expr):
        for c in expr.children:
            walk(c)
        # see if this is a magnitude cut
        # which is a lower limit on flux
        # so match for the pattern
        # Node_with_band_attribute > Literal
        if not isinstance(expr, Expr):
            return
        if expr.operator != '>':
            return
        if not isinstance(expr.children[1], Literal):
            return
        if not hasattr(expr.children[0], 'band'):
            return
        result.append(expr.children[0].band)
    walk(expr)
    return list(set(result))

# now we try to import a local version the file
def _local():

    import os
    import warnings

    local = os.path.join(os.path.dirname(__file__), 'local-%s' % os.path.basename(__file__))
    if local.endswith('.pyc'): local = local[:-1]

    if os.path.exists(local):
        warnings.warn(
             "Using local-target-selection.py is deprecated. \n" +
             ("Rename %s to my-targets.py, and add it to the commandline via " % local) +
             "--extra-target-definitions=my-targets.py" , stacklevel=2) 
        load(local)

def _export(g):
    global __all__

    # This will filter out names that does not appear to be a target type.
    import imaginglss.model.columnnames as columnnames
    blacklist = dir(columnnames)
    blacklist.extend(['load', 'WFLUX', 'GRZFLUX', 'SNRW1', 'SNRW2'])
    g1 = globals()
    for k in g.keys():
        if k in blacklist:
            continue
        if k in ['Min', 'Max']:
            continue
        if k.startswith('_'):
            continue
        g1[k] = g[k]
        item = g[k]
        item.name = k
        item.color_bands = _gather_color_bands(item)
        item.mag_bands = _gather_magnitude_bands(item)

        __all__.append(k)


_export(globals())
_local()
