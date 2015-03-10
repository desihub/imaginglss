# most code taken from https://raw.githubusercontent.com/dstndstn/tractor/master/projects/desi/common.py
# switch to astropy.wcs 
import fitsio
import os
from astropy import wcs
import numpy as np

# to run need to 
# export DUST_DIR=/project/projectdirs/desi/software/edison/dust/v0_0/

class anwcs_t(object):
    def __init__(self, filename, hduid):
        header = fitsio.read_header(filename, hduid)
        self.obj = wcs.WCS(header)
    def radec2pixelxy(self, ra, dec):
        input = np.array([ra, dec]).T
        x,y = self.obj.all_world2pix(input, 1).T
        return np.zeros_like(x), x, y
def radectolb(ra, dec):
    l, b = euler(ra, dec, 1)
    return l, b 
class SFDMap(object):
    # These come from Schlafly & Finkbeiner, arxiv 1012.4804v2, Table 6, Rv=3.1
    # but updated (and adding DES u) via email from Schlafly,
    # decam-data thread from 11/13/2014, "New recommended SFD coefficients for DECam."
    #
    # The coefficients for the four WISE filters are derived from Fitzpatrick 1999,
    # as recommended by Schafly & Finkbeiner, considered better than either the
    # Cardelli et al 1989 curves or the newer Fitzpatrick & Massa 2009 NIR curve
    # not vetted beyond 2 micron).
    # These coefficients are A / E(B-V) = 0.184, 0.113, 0.0241, 0.00910. 
    #
    extinctions = {
        'SDSS u': 4.239,
        'DES u': 3.995,
        'DES g': 3.214,
        'DES r': 2.165,
        'DES i': 1.592,
        'DES z': 1.211,
        'DES Y': 1.064,
        'WISE W1': 0.184,
        'WISE W2': 0.113,
        'WISE W3': 0.0241,
        'WISE W4': 0.00910,
        }

    def __init__(self, ngp_filename=None, sgp_filename=None, dustdir=None):
        if dustdir is None:
            dustdir = os.environ.get('DUST_DIR', None)
        if dustdir is not None:
            dustdir = os.path.join(dustdir, 'maps')
        else:
            dustdir = '.'
            print 'Warning: $DUST_DIR not set; looking for SFD maps in current directory.'
        if ngp_filename is None:
            ngp_filename = os.path.join(dustdir, 'SFD_dust_4096_ngp.fits')
        if sgp_filename is None:
            sgp_filename = os.path.join(dustdir, 'SFD_dust_4096_sgp.fits')
        if not os.path.exists(ngp_filename):
            raise RuntimeError('Error: SFD map does not exist: %s' % ngp_filename)
        if not os.path.exists(sgp_filename):
            raise RuntimeError('Error: SFD map does not exist: %s' % sgp_filename)
        self.north = fitsio.read(ngp_filename)
        self.south = fitsio.read(sgp_filename)
        self.northwcs = anwcs_t(ngp_filename, 0)
        self.southwcs = anwcs_t(sgp_filename, 0)

    @staticmethod
    def bilinear_interp_nonzero(image, x, y):
        H,W = image.shape
        x0 = np.floor(x).astype(int)
        y0 = np.floor(y).astype(int)
        # Bilinear interpolate, but not outside the bounds (where ebv=0)
        fx = np.clip(x - x0, 0., 1.)
        ebvA = image[y0,x0]
        ebvB = image[y0, np.clip(x0+1, 0, W-1)]
        ebv1 = (1.-fx) * ebvA + fx * ebvB
        ebv1[ebvA == 0] = ebvB[ebvA == 0]
        ebv1[ebvB == 0] = ebvA[ebvB == 0]

        ebvA = image[np.clip(y0+1, 0, H-1), x0]
        ebvB = image[np.clip(y0+1, 0, H-1), np.clip(x0+1, 0, W-1)]
        ebv2 = (1.-fx) * ebvA + fx * ebvB
        ebv2[ebvA == 0] = ebvB[ebvA == 0]
        ebv2[ebvB == 0] = ebvA[ebvB == 0]

        fy = np.clip(y - y0, 0., 1.)
        ebv = (1.-fy) * ebv1 + fy * ebv2
        ebv[ebv1 == 0] = ebv2[ebv1 == 0]
        ebv[ebv2 == 0] = ebv1[ebv2 == 0]
        return ebv

    def ebv(self, ra, dec):
        l,b = radectolb(ra, dec)
        ebv = np.zeros_like(l)
        N = (b >= 0)
        for wcs,image,cut in [(self.northwcs, self.north, N),
                              (self.southwcs, self.south, np.logical_not(N))]:
            # Our WCS routines are mis-named... the SFD WCSes convert 
            #   X,Y <-> L,B.
            if sum(cut) == 0:
                continue
            ok,x,y = wcs.radec2pixelxy(l[cut], b[cut])
            assert(np.all(ok == 0))
            H,W = image.shape
            assert(np.all(x >= 0.5))
            assert(np.all(x <= (W+0.5)))
            assert(np.all(y >= 0.5))
            assert(np.all(y <= (H+0.5)))
            ebv[cut] = SFDMap.bilinear_interp_nonzero(image, x-1., y-1.)
        return ebv

    def extinction(self, filts, ra, dec, get_ebv=False):
        ebv = self.ebv(ra, dec)
        factors = np.array([SFDMap.extinctions[f] for f in filts])
        rtn = factors[np.newaxis,:] * ebv[:,np.newaxis]
        if get_ebv:
            return ebv,rtn
        return rtn

from numpy import *

def euler(ai, bi, select=1, fk4=False):
   """
    NAME:
        EULER
    PURPOSE:
        Transform between Galactic, celestial, and ecliptic coordinates.
    EXPLANATION:
        Use the procedure ASTRO to use this routine interactively
   
    CALLING SEQUENCE:
         EULER, AI, BI, AO, BO, [ SELECT, /FK4, SELECT = ]
   
    INPUTS:
          AI - Input Longitude in DEGREES, scalar or vector.  If only two
                  parameters are supplied, then  AI and BI will be modified to
                  contain the output longitude and latitude.
          BI - Input Latitude in DEGREES
   
    OPTIONAL INPUT:
          SELECT - Integer (1-6) specifying type of coordinate transformation.
   
         SELECT   From          To        |   SELECT      From            To
          1     RA-Dec (2000)  Galactic   |     4       Ecliptic      RA-Dec
          2     Galactic       RA-DEC     |     5       Ecliptic      Galactic
          3     RA-Dec         Ecliptic   |     6       Galactic      Ecliptic
   
         If not supplied as a parameter or keyword, then EULER will prompt for
         the value of SELECT
         Celestial coordinates (RA, Dec) should be given in equinox J2000
         unless the /FK4 keyword is set.
    OUTPUTS:
          AO - Output Longitude in DEGREES
          BO - Output Latitude in DEGREES
   
    INPUT KEYWORD:
          /FK4 - If this keyword is set and non-zero, then input and output
                celestial and ecliptic coordinates should be given in equinox
                B1950.
          /SELECT  - The coordinate conversion integer (1-6) may alternatively be
                 specified as a keyword
    NOTES:
          EULER was changed in December 1998 to use J2000 coordinates as the
          default, ** and may be incompatible with earlier versions***.
    REVISION HISTORY:
          Written W. Landsman,  February 1987
          Adapted from Fortran by Daryl Yentis NRL
          Converted to IDL V5.0   W. Landsman   September 1997
          Made J2000 the default, added /FK4 keyword  W. Landsman December 1998
          Add option to specify SELECT as a keyword W. Landsman March 2003
   """

   n_params = 5
   select1 = select
   
   # ON_ERROR, 2
   
#   print 'Syntax - EULER, AI, BI, A0, B0, [ SELECT, /FK4, SELECT= ]'
#   print '    AI,BI - Input longitude,latitude in degrees'
#   print '    AO,BO - Output longitude, latitude in degrees'
#   print '    SELECT - Scalar (1-6) specifying transformation type'
   
   twopi = 2.0e0 * pi
   fourpi = 4.0e0 * pi
   deg_to_rad = 180.0e0 / pi
   
   #   J2000 coordinate conversions are based on the following constants
   #   (see the Hipparcos explanatory supplement).
   #  eps = 23.4392911111d              Obliquity of the ecliptic
   #  alphaG = 192.85948d               Right Ascension of Galactic North Pole
   #  deltaG = 27.12825d                Declination of Galactic North Pole
   #  lomega = 32.93192d                Galactic longitude of celestial equator
   #  alphaE = 180.02322d              Ecliptic longitude of Galactic North Pole
   #  deltaE = 29.811438523d            Ecliptic latitude of Galactic North Pole
   #  Eomega  = 6.3839743d              Galactic longitude of ecliptic equator
   
   if fk4:   
      equinox = '(B1950)'
      psi = array ([0.57595865315e0, 4.9261918136e0, 0.00000000000e0, 0.0000000000e0, 0.11129056012e0, 4.7005372834e0])
      stheta = array ([0.88781538514e0, -0.88781538514e0, 0.39788119938e0, -0.39788119938e0, 0.86766174755e0, -0.86766174755e0])
      ctheta = array([0.46019978478e0, 0.46019978478e0, 0.91743694670e0, 0.91743694670e0, 0.49715499774e0, 0.49715499774e0])
      phi = array([4.9261918136e0, 0.57595865315e0, 0.0000000000e0, 0.00000000000e0, 4.7005372834e0, 0.11129056012e0])
   else:   
      equinox = '(J2000)'
      psi = array([0.57477043300e0, 4.9368292465e0, 0.00000000000e0, 0.0000000000e0, 0.11142137093e0, 4.71279419371e0])
      stheta = array([0.88998808748e0, -0.88998808748e0, 0.39777715593e0, -0.39777715593e0, 0.86766622025e0, -0.86766622025e0])
      ctheta = array([0.45598377618e0, 0.45598377618e0, 0.91748206207e0, 0.91748206207e0, 0.49714719172e0, 0.49714719172e0])
      phi = array([4.9368292465e0, 0.57477043300e0, 0.0000000000e0, 0.00000000000e0, 4.71279419371e0, 0.11142137093e0])
      
   i = select - 1                         # IDL offset
   a = ai / deg_to_rad - phi[i]
   b = bi / deg_to_rad
   sb = sin(b) ;        cb = cos(b)
   cbsa = cb * sin(a)
   b = -stheta[i] * cbsa + ctheta[i] * sb
   bo = arcsin(minimum(b, 1.0e0)) * deg_to_rad

   a = arctan2(ctheta[i] * cbsa + stheta[i] * sb, cb * cos(a))
   ao = ((a + psi[i] + fourpi) % twopi) * deg_to_rad

   return (ao,bo)

if __name__ == '__main__':
    from model.datarelease import DataRelease
    dr = DataRelease()
    RA = dr.catalogue['RA']
    DEC = dr.catalogue['DEC']
    EXT = dr.catalogue['DECAM_EXTINCTION']
    m = SFDMap()
    EXT2 = m.extinction(['DES %s' % i for i in 'ugrizY'],
            RA, DEC)
    print EXT - EXT2
