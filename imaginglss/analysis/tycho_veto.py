"""
    veto objects based on Tycho catalogue. This
    is based on the email discussion at:

    Date: June 18, 2015 at 3:44:09 PM PDT
    To: decam-data@desi.lbl.gov
    Subject: decam-data Digest, Vol 12, Issue 29

    Current usage is to first create a Tycho catalogue
    object, then use DECAM_XXX functions to create a
    mask :

    >>> tycho = Tycho("pathto_tycho.fit")
    >>> mask = DECAM_LRG(tycho, coord)
    >>>

    The return values are True for 'preserve',
    False for 'reject'
"""

def BOSS_DR9(tycho, coord):
    """ Returns True for 'keep' """
    # BOSS DR9-11
    d, bmag, vmag = tycho.nearest(coord)
    b = bmag.clip(6, 11.5)
    R = (0.0802 * b ** 2 - 1.86 * b + 11.625) / 60. # 
    return d > R

BOSS_DR9.BIT = 1 << 0

def DECAM_LRG(tycho, coord):
    ra, dec = coord
    d, bmag, vmag = tycho.nearest(coord)
    R = 10 ** (3.5 - 0.15 * vmag) / 3600. 
    return d > R

DECAM_LRG.BIT = 1 << 1

DECAM_ELG = DECAM_LRG

DECAM_ELG.BIT = 1 << 2

def DECAM_QSO(tycho, coord):
    # I recommend not applying a bright star mask 
    ra, dec = coord
    return ra == ra

DECAM_QSO.BIT = 1 << 3

def DECAM_BGS(tycho, coord):
    ra, dec = coord
    d, bmag, vmag = tycho.nearest(coord)
    R = 10 ** (2.2 - 0.15 * vmag) / 3600. 
    return d > R

DECAM_BGS.BIT = 1 << 4
