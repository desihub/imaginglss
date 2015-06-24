

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

DECAM_INVVAR = C('DECAM_INVVAR').T
DECAM_MW_TRANSMISSION = C('DECAM_MW_TRANSMISSION').T

# 1-sigma limits
G_LIMIT = DECAM_INVVAR[1] ** -0.5 / DECAM_MW_TRANSMISSION[1]
Z_LIMIT = DECAM_INVVAR[4] ** -0.5 / DECAM_MW_TRANSMISSION[4]
R_LIMIT = DECAM_INVVAR[2] ** -0.5 / DECAM_MW_TRANSMISSION[2]

def LRG(sigma_r, sigma_z):
    LRG  = sigma_r * R_LIMIT < 10**((22.5-23.00    )/2.5)
    LRG &= sigma_z * Z_LIMIT < 10**((22.5-20.56    )/2.5)
    LRG &= sigma_z * Z_LIMIT < 10**((22.5-23.00+1.6)/2.5)
    return LRG

def ELG(sigma_g, sigma_r, sigma_z):
    ELG =  sigma_g * G_LIMIT < 10**((22.5-23.4-1.5+0.2)/2.5)
    ELG &= sigma_r * R_LIMIT < 10**((22.5-23.4        )/2.5)
    ELG &= sigma_z * Z_LIMIT < 10**((22.5-23.4+0.3    )/2.5)
    return ELG

def QSO(sigma_r, sigma_g):
    QSO &= sigma_r * R_LIMIT < 10**((22.5-23.00    )/2.5)
    QSO &= sigma_g * G_LIMIT < 10**((22.5-23.00-1.0)/2.5)
    return QSO

def BGS(sigma_r):
    BGS = sigma_r * R_LIMIT < 10**((22.5-19.5)/2.5)
    return BGS
