# Python routines to apply world-coordinate system transformations
# based on the tangent-plane projection.
# This code is designed to be light-weight, and specialized.
#
# The tangent-plane projection is described in
# http://fits.gsfc.nasa.gov/registry/tpvwcs/tpv.html
# and in:
# "Representations of celestial coordinates in FITS"
# -- Calabretta, M. R., and Greisen, E. W.,
# -- Astronomy & Astrophysics, 395, 1077-1122, 2002.
#

__author__ = "Yu Feng and Martin White"
__version__ = "0.9"
__email__  = "yfeng1@berkeley.edu or mwhite@berkeley.edu"


import numpy as N

def ang2pix(hdr,ra,dc):
    """
    ang2pix(hdr,ra,dc):
    Given a WCS header hdr (a Python dict) and numpy arrays of
    RA and DEC (in decimal degrees) returns the pixel numbers
    associated with the points.
    If a given RA,DEC lies outside the range, -1 is returned.
    """
    # Do a check that we are looking at a tangent plane projection,
    # else throw an exception.
    xx = N.zeros_like(ra)
    yy = N.zeros_like(dc)
    return( (xx,yy) )
    #


def pix2ang(hdr,xx,yy):
    """
    pix2ang(hdr,xx,yy):
    Given a WCS header hdr (a Python dict) and numpy arrays of
    pixel numbers (x and y) returns RA and DEC (in decimal degrees).
    If a given pixel lies outside the range, -1 is returned.
    """
    # Do a check that we are looking at a tangent plane projection,
    # else throw an exception.
    ra = N.zeros_like(xx)
    dc = N.zeros_like(yy)
    return( (ra,dc) )
    #
