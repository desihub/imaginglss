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

import numpy 
import math
__all__ = ['ang2pix', 'pix2ang']
def ang2pix(coord, CD, CRPIX, CRVAL):
    """
    TAN Transform from coord = (ra, dec) to pixel xy coordinate 
    according the the WCS header. Look up Section 5.1.3 of 
    http://www.aanda.org/articles/aa/pdf/2002/45/aah3860.pdf 
    
    Although the source code in 

        https://code.google.com/p/esutil/source/browse/trunk/esutil/wcsutil.py

    maybe a better explanation of what is done.

    coord = (ra, dec), RA and DEC (in decimal degrees, vectorized) 

    returns the pixel xy = (x, y). 

    Depending on the value of CRPIX, it can be 
        starting from 1 if the CRPIX1 and CRPIX2 in FITS header are used.
        starting from 0 if CRPIX1 - 1 and CRPIX2 - 1 are used (useful in numpy)

    CD is the tranformation matrix in CD1_1, CD1_2, CD2_1, CD_2_2,
    CRPIX is the center of pixels as of original FITS header (starting from 1), CRPIX1, CRPIX2
    CRVAL is the center of RA/DEC as of original FITS header. CRVAL1, CRVAL2
    Obviously PV distortion is not supported.

    No checking is performed if a given RA, DEC lies outside the range.
    """

    coord = numpy.array(coord, dtype='f8').copy()
    xy = numpy.empty_like(coord)
    # watch out, this may be wrong if the matrix is not diagonal
    matrix = numpy.linalg.inv(numpy.array(CD).reshape(2, 2))

    r = CreateRotationMatrix(CRVAL[0], CRVAL[1])
    ra, dec = Rotate(r, coord[0], coord[1]) 
    ra *= numpy.pi / 180.
    dec *= numpy.pi / 180.

    rdiv = 180. / numpy.pi / numpy.tan(dec)
    xy[0] = rdiv * numpy.sin(ra)
    xy[1] = -rdiv * numpy.cos(ra)

    xy = matrix.dot(xy)
    xy += numpy.array(CRPIX).reshape(2, 1)

    return  xy


def pix2ang(xy, CD, CRPIX, CRVAL):
    """
        This is the invert transformation of ang2pix. (TAN)

        xy = (x, y): coordinate in pixels 

        Depending on the value of CRPIX, the offset of xy can be 
            starting from 1 if the CRPIX1 and CRPIX2 in FITS header are used.
            starting from 0 if CRPIX1 - 1 and CRPIX2 - 1 are used (useful in numpy)

        Returns coord = (RA, DEC), in degrees.
    """

    xy = numpy.array(xy, dtype='f8').copy()
    coord = numpy.empty_like(xy)
    # watch out, this may be wrong if the matrix is not diagonal
    matrix = numpy.array(CD).reshape(2, 2)

    xy -= numpy.array(CRPIX).reshape(2, 1)
    xy = matrix.dot(xy)

    rinv = numpy.einsum('ij,ij->j', xy, xy) ** -0.5
    rinv *= 180.0 / numpy.pi
    coord[1] = numpy.arctan(rinv)

    coord[0] = numpy.arctan2(xy[0], -xy[1])

    r = CreateRotationMatrix(CRVAL[0], CRVAL[1])
    coord[:] = Rotate(r.T, coord[0], coord[1]) 
    coord *= 180 / numpy.pi

    return coord


##############################
#  Private routines, taken from esutils.py
#
def CreateRotationMatrix(native_longpole, native_latpole):
    """ create rotation matrix for tan. adapted from esutil.py """
    longpole = 180.
    d2r = numpy.pi / 180.
    # If Theta0 = 90 then CRVAL gives the coordinates of the origin in the
    # native system.   This must be converted (using Eq. 7 in Greisen &
    # Calabretta with theta0 = 0) to give the coordinates of the North
    # pole (longitude_p, latitude_p)

    # Longpole is the longitude in the native system of the North Pole in
    # the standard system (default = 180 degrees).
    sp = math.sin(longpole*d2r)
    cp = math.cos(longpole*d2r)

    sa = math.sin(native_longpole * d2r)
    ca = math.cos(native_longpole * d2r)
    sd = math.sin(native_latpole * d2r)
    cd = math.cos(native_latpole * d2r)

    # calculate rotation matrix

    # esutils has this transpolsed
    r = numpy.array([[-sa*sp - ca*cp*sd,   sa*cp - ca*sp*sd, ca*cd ] ,
                     [ ca*sp - sa*cp*sd , -ca*cp - sa*sp*sd, sa*cd ] ,
                     [ cp*cd           ,   sp*cd           , sd    ] ],
                    dtype='f8')

    # we transpose it back 
    return r.T.copy()

def Rotate(r, longitude, latitude):
    """
    Apply a rotation matrix to the input longitude and latitude
    inputs must be numpy arrays; stolen from esutils.py
    """
    latitude = latitude * (numpy.pi / 180)
    longitude = longitude * (numpy.pi / 180)
    l = numpy.cos(latitude)*numpy.cos(longitude)
    m = numpy.cos(latitude)*numpy.sin(longitude)
    n = numpy.sin(latitude)

    lmn = numpy.array([l, m, n])

    # find solution to the system of equations and put it in b
    # Can't use matrix notation in case l,m,n are arrays

    lmn2 = r.dot(lmn)

    b0, b1, b2 = lmn2

    # Account for possible roundoff
    numpy.clip(b2, -1, 1, b2)

    lat_new = numpy.arcsin(b2)* (180. / numpy.pi)
    lon_new = numpy.arctan2(b1, b0) * (180. / numpy.pi)

    return lon_new, lat_new


if __name__ == '__main__':
    # perform some tests
    def compare(ra0, dec0, ra, dec):
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

        print 'transforming', ra, dec, 'at', ra0, dec0
        print 'roundtrip', ra, dec
        print 'astropy has', astropy
        print 'we have    ', ours.T
        return ours - astropy

    def test():
        compare(30, 30, 31., 30)
        compare(30, 30, 30., 30)
        compare(30, 30, 30., 31)
        compare(-180, 30, -181., 31)
        compare(0., 80, 0., 80.9)
        # we still have an issue at the north pole 
        # is this the right name for DEC=90? -YF
        compare(0., 90, 0., 89.9)
    test()
