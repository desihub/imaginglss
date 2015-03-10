# Python routines to apply world-coordinate system transformations
# based on the simpler ZEA projection used SFD98 milky way dust map
# This code is designed to be light-weight, and specialized.
#
# We handle two special cases 
# CRVAL = -90 -90 (south) and 90 90 (north)
# The ZEA projection is described in
# http://fits.gsfc.nasa.gov/registry/tpvwcs/tpv.html
# and in:
# "Representations of celestial coordinates in FITS"
# -- Calabretta, M. R., and Greisen, E. W.,
# -- Astronomy & Astrophysics, 395, 1077-1122, 2002.
# and the source code in 
#   https://code.google.com/p/esutil/source/browse/trunk/esutil/wcsutil.py
# is also a very useful reference.
#
"""
example header:
SIMPLE  =                    T / Written by IDL:  Sun May  9 11:19:26 1999
BITPIX  =                  -32 /
NAXIS   =                    2 /
NAXIS1  =                 4096 /
NAXIS2  =                 4096 /
DATE    = '1999-05-09'         / Creation date (CCYY-MM-DD) of FITS header
OBJECT  = 'E(B-V)  '           /
BUNIT   = 'mag     '           /
CRPIX1  =              2048.50 / 1-indexed X pixel number of pole
CRVAL1  =              90.0000 /
CTYPE1  = 'GLON-ZEA'           / X=sqrt(1-NSGP*sin(b))*cos(l)*SCALE
CRPIX2  =              2048.50 / 1-indexed Y pixel number of pole
CRVAL2  =              90.0000 /
CTYPE2  = 'GLAT-ZEA'           / Y=-NSGP*sqrt(1-NSGP*sin(b))*sin(l)*SCALE
CD1_1   =     -0.0395646818624 /
CD1_2   =              0.00000 /
CD2_1   =              0.00000 /
CD2_2   =      0.0395646818624 /
LONPOLE =                  180 /
LAM_NSGP=                    1 / NSGP=+1 for north polar, =-1 for south polar
LAM_SCAL=                 2048 / SCALE=number of pixels from b=0 to b=90 deg
AUTHOR  = 'David J. Schlegel, Douglas P. Finkbeiner, and Marc Davis' /
"""
from __future__ import print_function

__author__ = "Yu Feng and Martin White"
__version__ = "0.9"
__email__  = "yfeng1@berkeley.edu or mjwhite@lbl.gov"
__all__ = ['ang2pix','pix2ang','ang2pix_hdr','pix2ang_hdr']



import numpy 
import math

def ang2pix_hdr(coord,hdr,zero_offset=True):
    """
    A convenience function for calling ang2pix, with the transformations
    described in a dictionary in the usual FITS header style for WCS
    transformations.  Upon input coord = (ra, dec), RA and DEC (in decimal
    degrees, vectorized) and the routine returns pixel numbers.
    If zero_offset is True, the routine returns 0-indexed pixel coordinates
    (useful in Python or C) while if it is False pixels run from 1 (as in
    Fortran).
    """
    if ('ZEA' not in hdr['CTYPE1'])|('ZEA' not in hdr['CTYPE2']):
        raise RuntimeError,"Not a zea projection."
    cd, crpix, crval = parse_header(hdr, zero_offset)
    return(ang2pix(coord,cd,crpix,crval))
    #

def pix2ang_hdr(xy,hdr,zero_offset=True):
    """
    A convenience function for calling pix2ang, with the transformations
    described in a dictionary in the usual FITS header style for WCS
    transformations.  Given input pixel numbers xy the routine returns
    coordinates (RA,DEC) in decimal degrees.
    If zero_offset is True, the routine takes 0-indexed pixel coordinates
    (useful in Python or C) while if it is False pixels run from 1 (as in
    Fortran).
    """
    if ('ZEA' not in hdr['CTYPE1'])|('ZEA' not in hdr['CTYPE2']):
        raise RuntimeError,"Not a zea plane projection."
    cd, crpix, crval = parse_header(hdr, zero_offset)
    return(pix2ang(xy,cd,crpix,crval))
    #

def parse_header(hdr, zero_offset):
    """
    If zero_offset is True, the routine takes 0-indexed pixel coordinates
    (useful in Python or C) while if it is False pixels run from 1 (as in
    """
    # Check to see whether the "hdr" dictionary contains the necessary
    # keywords.
    if ('CTYPE1' not in hdr)|('CTYPE2' not in hdr)|\
       ('CRVAL1' not in hdr)|('CRVAL2' not in hdr)|\
       ('CRPIX1' not in hdr)|('CRPIX2' not in hdr)|\
       ('LAM_NSGP' not in hdr)|('LAM_SCALE' not in hdr):
        raise RuntimeError,"Unable to parse header."
    crpix = numpy.array([hdr['CRPIX1'],hdr['CRPIX2']])
    nsgp = int(hdr['LAM_NSGP'])
    scale = int(hdr['LAM_SCALE'])

    if zero_offset:
        crpix -= 1
    return scale, crpix, nsgp
               


def ang2pix(coord,SCALE,CRPIX,NSGP):
    """
    ZEA Transform from coord = (ra, dec) to pixel xy coordinate 
    according the the WCS header. Look up Section 5.?.? of 
    http://www.aanda.org/articles/aa/pdf/2002/45/aah3860.pdf 
    Although the source code in 
      https://code.google.com/p/astropysics/source/browse/astropysics/extinct.py?r=93bcf1e49124f3fe06f8369ac290fe7d8a8f80fc

    maybe a better explanation of what is done.
    coord = (ra, dec), RA and DEC (in decimal degrees, vectorized) 
    returns the pixel xy = (x, y). 

    CD is the tranformation matrix in CD1_1, CD1_2, CD2_1, CD_2_2,

    RA/DEC at CRVAL shall map exactly to x/y at CRPIX.

    Obviously PV distortion is not supported.
    No checking is performed if a given RA, DEC lies outside the range.
    """
    coord = numpy.array(coord, dtype='f8').copy()
    xy    = numpy.empty_like(coord)
    # watch out, this may be wrong if the matrix is not diagonal

    l, b = coord * (numpy.pi / 180.)
       
    #project from galactic longitude/latitude to lambert pixels (see SFD98)
    x = numpy.cos(l) * (1 - NSGP * numpy.sin(b))**0.5
    y = - NSGP *numpy.sin(l) * (1 - NSGP * numpy.sin(b))**0.5
    #now remap indecies - numpy arrays have y and x convention switched
    xy[0] = x
    xy[1] = y

    xy *= SCALE
    xy += numpy.array(CRPIX).reshape(2, 1)
    return(xy)
    #


def pix2ang(xy,CD,CRPIX,CRVAL):
    """
    This is the invert transformation of ang2pix. (TAN)
    xy = (x, y): coordinate in pixels
    Returns coord = (RA, DEC), in decimal degrees.

    RA/DEC at CRVAL shall map exactly to x/y at CRPIX.
    """
    xy    = numpy.array(xy, dtype='f8').copy()
    coord = numpy.empty_like(xy)
    # watch out, this may be wrong if the matrix is not diagonal
    matrix = numpy.array(CD).reshape(2, 2)

    xy  -= numpy.array(CRPIX).reshape(2, 1)
    xy   = matrix.dot(xy)
    rinv = numpy.einsum('ij,ij->j',xy,xy)
    # this will give a reasonable coord[1] at pole
    rinv.clip(1e-28, out=rinv)
    rinv **= -0.5
    rinv  *= 180.0 / numpy.pi

    coord[1] = numpy.arctan(rinv)
    coord[0] = numpy.arctan2(xy[0],-xy[1])
    coord   *= 180 / numpy.pi
    coord[0]%= 360.
    return(coord)

if __name__ == '__main__':
    # perform some tests
    def compare(NSGP, ra, dec):
        from astropy import wcs
        header = dict(
            CTYPE1  = 'GLON-ZEA',#           / ZEA 
            CTYPE2  = 'GLAT-ZEA', #           / ZEA
            CRPIX1  =               2048.5, # / Reference x 
            CRPIX2  =               2048.5, # / Reference y 
            CD1_1   = NSGP*-0.0395646818624, # / CD matrix
            CD1_2   =                   0., # / CD matrix
            CD2_1   =                   0., # / CD matrix
            CD2_2   = NSGP*0.0395646818624, # / CD matrix
            LONPOLE = 180,
        )
        header['CRVAL1']  =     NSGP * 90. # / Reference RA                                   
        header['CRVAL2']  =     NSGP * 90.# / Reference Dec                                  
        q = wcs.WCS(header)
        
        ra = [ra]
        dec = [dec]
        astropy = q.all_world2pix(numpy.array((ra, dec)).T, 1)
        ours = ang2pix((ra, dec), 
                SCALE=2048,
                CRPIX=(2048.5,2048.5),
                NSGP=NSGP)
        back = pix2ang(ours, 
                CD=(0.0395646818624,0., 0.,-0.0395646818624,),
                CRPIX=(1800.5,1800.5),
                CRVAL=(ra, dec))

        print('transforming', ra, dec, 'at', ra, dec)
        print('roundtrip', back.T)
        print('astropy has', astropy)
        print('we have    ', ours.T)
        return ours - astropy

    def test():
        compare(1, 31., -30)
        compare(1, 30., -30)
        compare(1, 30., -31)
        compare(1, -181., -31)
    test()
