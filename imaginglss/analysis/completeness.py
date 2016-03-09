

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

    Put your own complete area definitions in

        local-completeness.py 

    of the same path.

Remember to append the name of the object type to __all__ variable
with __all__.append("ObjectType")

"""
from imaginglss.utils.npyquery import Column as C
__all__ = []

DECAM_DEPTH = C('DECAM_DEPTH')
DECAM_MW_TRANSMISSION = C('DECAM_MW_TRANSMISSION')

# 1-sigma limits
G_LIMIT = DECAM_DEPTH[1] ** -0.5 / DECAM_MW_TRANSMISSION[1]
Z_LIMIT = DECAM_DEPTH[4] ** -0.5 / DECAM_MW_TRANSMISSION[4]
R_LIMIT = DECAM_DEPTH[2] ** -0.5 / DECAM_MW_TRANSMISSION[2]

def LRG(sigma):
    """ Create the completeness cut for LRG.

        Parameters
        ----------
        sigma : dict 
            'r', 'z', 'g' members specifies the confidents in that
            band. ( :math:`n \\sigma`)

    """
    LRG  = sigma['r'] * R_LIMIT < 10**((22.5-23.00    )/2.5)
    LRG &= sigma['z'] * Z_LIMIT < 10**((22.5-20.56    )/2.5)
    LRG &= sigma['z'] * Z_LIMIT < 10**((22.5-23.00+1.6)/2.5)
    return LRG
LRG.bands = 'rz'
__all__.append("LRG")

def ELG(sigma):
    """ Create the completeness cut for ELG.

        Parameters
        ----------
        sigma : dict 
            'r', 'z', 'g' members specifies the confidents in that
            band. ( :math:`n \\sigma`)

    """
    ELG = sigma['r'] * R_LIMIT < 10**((22.5-23.4        )/2.5)
    ELG &= sigma['z'] * Z_LIMIT < 10**((22.5-23.4+0.3    )/2.5)
    ELG &= sigma['g'] * G_LIMIT < 10**((22.5-23.4-1.5+0.2)/2.5)
    ELG &= sigma['g'] * G_LIMIT < 10**((22.5-23.4-1.2+0.3)/2.5)
    return ELG
ELG.bands = 'grz'
__all__.append("ELG")

def QSO(sigma):
    """ Create the completeness cut for QSO.

        Parameters
        ----------
        sigma : dict 
            'r', 'z', 'g' members specifies the confidents in that
            band. ( :math:`n \\sigma`)

    """
    QSO =  sigma['r'] * R_LIMIT < 10**((22.5-23.00    )/2.5)
    QSO &= sigma['g'] * G_LIMIT < 10**((22.5-23.00-1.0)/2.5)
    return QSO
QSO.bands = 'gr'
__all__.append("QSO")

def BGS(sigma):
    """ Create the completeness cut for BGS.

        Parameters
        ----------
        sigma : dict 
            'r', 'z', 'g' members specifies the confidents in that
            band. ( :math:`n \\sigma`)

    """
    BGS = sigma['r'] * R_LIMIT < 10**((22.5-19.5)/2.5)
    return BGS
BGS.bands = 'r'
__all__.append("BGS")

# now we try to import a local version the file
def _local():
    import os
    local = os.path.join(os.path.dirname(__file__), 'local-%s' % os.path.basename(__file__))
    if local.endswith('.pyc'): local = local[:-1]
    if os.path.exists(local):
        script = open(local, 'r').read()
        exec(compile(script, local, 'exec'), globals())
    del globals()['_local']
_local()
