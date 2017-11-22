"""
    veto objects based on a star catalogue.
    
    The "DESI_" vetos are based on Kitanidis++17 in prep.
    
    The "DECAM_" vetos are based on the email discussion at:

    Date: June 18, 2015 at 3:44:09 PM PDT
    To: decam-data@desi.lbl.gov
    Subject: decam-data Digest, Vol 12, Issue 29

    These objects takes a decals object and calculates the
    center and rejection radius for the catalogue in degrees.

    Note : The convention for veto flags is True for 'reject',
    False for 'preserve'.

    >>>
"""

def DESI_TYCHO(decals):
    tycho = decals.tycho
    vtmag = tycho['VTMAG']
    R = 10 ** (4.1 - 0.2 * vtmag) / 3600. 
    return tycho['RA'], tycho['DEC'], R

def DESI_WISE(decals):
    wise = decals.wise
    W1mpro = wise['W1MPRO']
    radius = 10 ** (3.0 - 0.1 * W1mpro) 
    radius[W1mpro < 2] = 1100.
    R = radius / 3600.
    return wise['RA'], wise['DEC'], R

def DECAM_LRG(decals):
    tycho = decals.tycho
    vtmag = tycho['VTMAG']
    R = 10 ** (3.5 - 0.15 * vtmag) / 3600. 
    return tycho['RA'], tycho['DEC'], R

DECAM_ELG = DECAM_LRG

def DECAM_QSO(decals):
    tycho = decals.tycho
    vtmag = tycho['VTMAG']
    # "I recommend not applying a bright star mask" -- D. Schlegal
    return tycho['RA'], tycho['DEC'], vtmag - vtmag

def DECAM_BGS(decals):
    tycho = decals.tycho
    vtmag = tycho['VTMAG']
#    R = 10 ** (2.2 - 0.15 * vtmag) / 3600. 
    # The above mask is not conservative enough. Ellie suggests:
    R = 10 ** (3.3 - 0.15 * vtmag) / 3600.
    return tycho['RA'], tycho['DEC'], R

def BOSS_DR9(decals):
    tycho = decals.tycho
    bmag = tycho['BMAG']
    # BOSS DR9-11
    b = bmag.clip(6, 11.5)
    R = (0.0802 * b ** 2 - 1.86 * b + 11.625) / 60. # 
    return tycho['RA'], tycho['DEC'], R

def EBOSS_V6(decals):
    wise = decals.wise
    W1mpro = wise['W1MPRO']
    # first do arcsecs then convert to degrees
    radius = 1397.5 - 569.34 * W1mpro + 79.88 * W1mpro ** 2 - 3.75 * W1mpro ** 3
    radius[W1mpro < 2] = 550.

    R = radius / 3600.
    return wise['RA'], wise['DEC'], R