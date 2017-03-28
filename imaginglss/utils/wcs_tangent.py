"""
Python routines to apply world-coordinate system transformations
based on the tangent-plane projection.
This code is designed to be light-weight, and specialized.

The tangent-plane projection is described in
http://fits.gsfc.nasa.gov/registry/tpvwcs/tpv.html
and in:
"Representations of celestial coordinates in FITS",
Calabretta, M. R., and Greisen, E. W.,
Astronomy & Astrophysics, 395, 1077-1122, 2002.

The source code in
   https://code.google.com/p/esutil/source/browse/trunk/esutil/wcsutil.py

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
    Apply forward (ra,dec)->(x,y) TAN transformation
    
    A convenience function for calling ang2pix, with the transformations
    described in a dictionary in the usual FITS header style for WCS
    transformations.  
    
    Parameters
    ----------
    coord  : array_like 
        (ra, dec), RA and DEC (in decimal degrees, vectorized) 
    hdr    :   dict
        WCS header as a dictionary
    zero_offset : boolean, optional,
        If True, the routine returns 0-indexed pixel coordinates
        (useful in Python or C) while if it is False pixels run from 1 (as in
        Fortran).

    Returns
    -------
    xy    : array_like
        xy = (x, y) pixel numbers.
    
    """
    if ('TAN' not in hdr['CTYPE1'])|('TAN' not in hdr['CTYPE2']):
        raise RuntimeError("Not a tangent plane projection.")
    cd, crpix, crval = parse_header(hdr, zero_offset)
    return(ang2pix(coord,cd,crpix,crval))
    #

def pix2ang_hdr(xy,hdr,zero_offset=True):
    """
    Apply backward (x,y)->(ra,dec) TAN transformation
    
    See :py:meth:`ang2pix_hdr`
    """
    if ('TAN' not in hdr['CTYPE1'])|('TAN' not in hdr['CTYPE2']):
        raise RuntimeError("Not a tangent plane projection.")
    cd, crpix, crval = parse_header(hdr, zero_offset)
    return(pix2ang(xy,cd,crpix,crval))
    #

def parse_header(hdr, zero_offset):
    """
    Parse a WCS header to arguments of pix2ang and ang2pix.

    Parameters
    ----------
    hdr         : dict
        WCS header
    zero_offset : boolean
        If zero_offset is True, the routine assumes 0-indexed pixel coordinates
        (useful in Python or C) while if it is False pixels run from 1 
        (as in Fortran and Julia)

    Returns
    -------
    cd    : array_like
        Transformation matrix (2, 2)
    crpix : array_like
        Centeral pixel number, compensated for `zero_offset`
    crval : array_like
        Centeral RA, DEC

    Raises
    ------
    RuntimeError:
        if the header does not contain enough fields.

    """
    # Check to see whether the "hdr" dictionary contains the necessary
    # keywords.
    if ('CTYPE1' not in hdr)|('CTYPE2' not in hdr)|\
       ('CRVAL1' not in hdr)|('CRVAL2' not in hdr)|\
       ('CRPIX1' not in hdr)|('CRPIX2' not in hdr)|\
       ('CD1_1'  not in hdr)|('CD1_2'  not in hdr)|\
       ('CD2_1'  not in hdr)|('CD2_2'  not in hdr):
        raise RuntimeError("Unable to parse header.")
    cd    = numpy.array([hdr['CD1_1'],hdr['CD1_2'],hdr['CD2_1'],hdr['CD2_2']])
    crpix = numpy.array([hdr['CRPIX1'],hdr['CRPIX2']])
    crval = numpy.array([hdr['CRVAL1'],hdr['CRVAL2']])
    if zero_offset:
        crpix -= 1
    return cd, crpix, crval
               


def ang2pix(coord,CD,CRPIX,CRVAL):
    """
    Convert RA, DEC to x,y, with TAN transformation

    Obviously PV distortion is not supported.

    No checking is performed if a given RA, DEC lies outside the range.

    Parameters
    ----------
    coord : array_like
        coord = (RA, DEC), RA and DEC (in decimal degrees, vectorized) 
    CD    : array_like
        transformation matrix (2, 2)
    CRPIX : array_like
        center pixel number of (x, y), compensated by offset.
    CRVAL : array_like
        center coordinate of (RA, DEC), in degrees.

    Notes
    -----
    Look up Section 5.?.? of 
    http://www.aanda.org/articles/aa/pdf/2002/45/aah3860.pdf 
    
    The code large follows implementation at
    https://code.google.com/p/esutil/source/browse/trunk/esutil/wcsutil.py
    """
    coord = numpy.array(coord, dtype='f8').copy()
    view = coord
    view = view.reshape(2, -1)
    CRVAL = numpy.array(CRVAL).reshape(2, -1)
    CRPIX = numpy.array(CRPIX).reshape(2, -1)

    xy    = numpy.empty_like(view)

    # watch out, this may be wrong if the matrix is not diagonal
    matrix = numpy.linalg.inv(numpy.array(CD).reshape(2, 2, -1).transpose((2, 0, 1)))

    ra, dec = native_transform(CRVAL[0], CRVAL[1], view[0], view[1], inverted=False)

    ra *= numpy.pi / 180.
    dec *= numpy.pi / 180.

    rdiv = 180. / numpy.pi / numpy.tan(dec)
    xy[0] = rdiv * numpy.sin(ra)
    xy[1] = -rdiv * numpy.cos(ra)

    xy = numpy.einsum('imn,ni->mi', matrix, xy)
    xy += numpy.array(CRPIX).reshape(2, -1)
    return xy.reshape(coord.shape)
    #
def native_transform(native_longpole, native_latpole, longitude, latitude, inverted=False):
    
    longpole = 180.
    d2r = numpy.pi / 180.
    # If Theta0 = 90 then CRVAL gives the coordinates of the origin in the
    # native system.   This must be converted (using Eq. 7 in Greisen &
    # Calabretta with theta0 = 0) to give the coordinates of the North
    # pole (longitude_p, latitude_p)

    # Longpole is the longitude in the native system of the North Pole in
    # the standard system (default = 180 degrees).
    sp = numpy.sin(longpole*d2r)
    cp = numpy.cos(longpole*d2r)

    sa = numpy.sin(native_longpole * d2r)
    ca = numpy.cos(native_longpole * d2r)
    sd = numpy.sin(native_latpole * d2r)
    cd = numpy.cos(native_latpole * d2r)

    # calculate rotation matrix

    # esutils has this transpolsed
    # we transpose it back 
    r = numpy.array([[-sa*sp - ca*cp*sd,   sa*cp - ca*sp*sd, ca*cd ] ,
                     [ ca*sp - sa*cp*sd , -ca*cp - sa*sp*sd, sa*cd ] ,
                     [ cp*cd           ,   sp*cd           , sd    ] ],
                    dtype='f8').transpose((2, 1, 0))

    
    latitude = latitude * (numpy.pi / 180)
    longitude = longitude * (numpy.pi / 180)
    x = numpy.cos(latitude)*numpy.cos(longitude)
    y = numpy.cos(latitude)*numpy.sin(longitude)
    z = numpy.sin(latitude)

    xyz = numpy.array([x, y, z])

    # find solution to the system of equations and put it in b
    # Can't use matrix notation in case l,m,n are arrays

    if not inverted:
        lmn2 = numpy.einsum('imn,ni->mi', r, xyz)
    else:
        lmn2 = numpy.einsum('inm,ni->mi', r, xyz)
    
    b0, b1, b2 = lmn2

    # Account for possible roundoff
    b2 = numpy.clip(b2, -1, 1)

    lat_new = numpy.arcsin(b2)* (180. / numpy.pi)
    lon_new = numpy.arctan2(b1, b0) * (180. / numpy.pi)

    return lon_new, lat_new

def pix2ang(xy,CD,CRPIX,CRVAL):
    """
    Convert x, y to RA, DEC with TAN transformation.

    See :py:meth:`ang2pix`
    
    """
    xy    = numpy.array(xy, dtype='f8').copy()
    coord = numpy.empty_like(xy)
    view = coord

    view = view.reshape(2, -1)
    xy = xy.reshape(2, -1)
    CRVAL = numpy.array(CRVAL).reshape(2, -1)
    CRPIX = numpy.array(CRPIX).reshape(2, -1)

    # watch out, this may be wrong if the matrix is not diagonal
    matrix = numpy.array(CD).reshape(2, 2, -1).transpose((2, 0, 1)).copy()

    xy  -= numpy.array(CRPIX).reshape(2, -1)
    xy   = numpy.einsum('imn,ni->mi', matrix, xy)
    rinv = numpy.einsum('mi,mi->i',xy,xy)

    # this will give a reasonable coord[1] at pole
    rinv.clip(1e-28, out=rinv)
    rinv **= -0.5
    rinv  *= 180.0 / numpy.pi

    view[1] = numpy.arctan(rinv)
    view[0] = numpy.arctan2(xy[0],-xy[1])
    view   *= 180 / numpy.pi
    ra, dec = native_transform(CRVAL[0], CRVAL[1], view[0], view[1], inverted=True)
    view[0]%= 360.
    return coord

if __name__ == '__main__':
    # perform some tests
    def compare(ra0, dec0, ra, dec):
        from numpy.testing import assert_allclose
        from astropy import wcs
        header = dict(
            CTYPE1  = 'RA---TAN',#           / TANgent plane
            CTYPE2  = 'DEC--TAN', #           / TANgent plane
            CRPIX1  =               1800.5, # / Reference x 
            CRPIX2  =               1800.5, # / Reference y 
            CD1_1   = -7.27777777777778E-05, # / CD matrix
            CD1_2   =                   0., # / CD matrix
            CD2_1   =                   0., # / CD matrix
            CD2_2   = 7.27777777777778E-05, # / CD matrix
        )
        header['CRVAL1']  =     ra0 # / Reference RA                                   
        header['CRVAL2']  =     dec0 # / Reference Dec                                  
        q = wcs.WCS(header)
        
        ra = [ra]
        dec = [dec]
        astropy = q.all_world2pix(numpy.array((ra, dec)).T, 1)
        ours = ang2pix((ra, dec), 
                CD=(-7.27777777777778E-05,0., 0.,7.27777777777778E-05,),
                CRPIX=(1800.5,1800.5),
                CRVAL=(ra0, dec0))
        back = pix2ang(ours, 
                CD=(-7.27777777777778E-05,0., 0.,7.27777777777778E-05,),
                CRPIX=(1800.5,1800.5),
                CRVAL=(ra0, dec0))

        print('transforming', ra, dec, 'at', ra0, dec0)
        print('roundtrip', back.T)
        print('astropy has', astropy)
        print('we have    ', ours.T)
        assert_allclose(ours.T, astropy, rtol=1e-9)

    def test():
        # assert vector input runs.
        ra = [30] * 2
        dec = [30] * 2
        ra0 = [40] * 2
        dec0 = [50] * 2

        ours = ang2pix((ra, dec), 
                CD=(-7.27777777777778E-05,0., 0.,7.27777777777778E-05,),
                CRPIX=(1800.5,1800.5),
                CRVAL=(ra0, dec0))
        back = pix2ang(ours, 
                CD=(-7.27777777777778E-05,0., 0.,7.27777777777778E-05,),
                CRPIX=(1800.5,1800.5),
                CRVAL=(ra0, dec0))

        compare(30, 30, 31., 30)
        compare(30, 30, 30., 30)
        compare(30, 30, 30., 31)
        compare(-180, 30, -181., 31)
        compare(0., 80, 0., 80.9)
        # we still have an issue at the north pole 
        # is this the right name for DEC=90? -YF
        compare(0., 90, 0., 89.9)

    test()
