"""
    Target Selection

    https://desi.lbl.gov/trac/wiki/TargetSelection

    A collection of helpful (static) methods to check whether an object's
    flux passes a given selection criterion (e.g. LRG, ELG or QSO).


"""

__author__ = "Yu Feng and Martin White"
__version__ = "1.0"
__email__  = "yfeng1@berkeley.edu or mjwhite@lbl.gov"

# * is evil but exactly what we want to do here
from imaginglss.model.columnnames import *

__all__ = []

def load(path=None):
    """ Recursively load target definitions in the path."""
    import os
    if path is None:
        # load the files in the definitions directory
        path = os.path.join(os.path.dirname(__file__), 'definitions')
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

def _is_target_definition(expr):
    """ an expr is a target definition, iif it has mag and color cuts
    """
    from imaginglss.utils.npyquery import Expr
    if not isinstance(expr, Expr): return False
    magbands = _gather_magnitude_bands(expr)
    colorbands = _gather_color_bands(expr)
    return len(magbands) > 0 and len(colorbands) > 0

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
    g1 = globals()
    for k in g.keys():
        if k.startswith('_'):
            continue
        item = g[k]
        if not _is_target_definition(item):
            continue
        item.name = k
        item.color_bands = _gather_color_bands(item)
        item.mag_bands = _gather_magnitude_bands(item)

        # register to global, and export in __all__
        g1[k] = item
        __all__.append(k)

load(None)

# backward compatable
_local()
