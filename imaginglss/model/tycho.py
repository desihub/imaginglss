import numpy
from ..utils import fits
from scipy.spatial import cKDTree

__author__ = "Yu Feng and Martin White"
__version__ = "1.0"
__email__  = "yfeng1@berkeley.edu or mjwhite@lbl.gov"

def radec2pos(ra, dec):
    """ converting ra dec to position on a unit sphere.
        ra, dec are in degrees.
    """
    pos = numpy.empty(len(ra), dtype=('f8', 3))
    ra = ra * (numpy.pi / 180)
    dec = dec * (numpy.pi / 180)
    pos[:, 2] = numpy.sin(dec)
    pos[:, 0] = numpy.cos(dec) * numpy.sin(ra)
    pos[:, 1] = numpy.cos(dec) * numpy.cos(ra)
    return pos

class Tycho(object):
    """ Representing a Tycho catalogue that
        is used to veto objects near stars 
    """
    def __init__(self, path):
        data = fits.read_table(path)
        self.bmag = data['BMAG'].copy()
        self.vmag = data['VMAG'].copy()
        self.varflag = data['VARFLAG'].copy()
        pos = radec2pos(data['RA'], data['DEC'])
        self.tree = cKDTree(pos)

    def nearest(self, coord):
        """
        Query Tycho catalogue, returns information about
        the nearest object. The distance is in degrees.
    
        Parameters
        ----------
        coord : array_like
            (ra, dec) in degrees

        Returns
        -------
        d, bmag, vmag: 
            the distance to each object in coord, and the
            mag in b and v.
        """
        pos = radec2pos(coord[0], coord[1])
        d, i = self.tree.query(pos, k=1)
        d = 2 * numpy.arcsin(d * 0.5) * 180 / numpy.pi
        return d, self.bmag[i], self.vmag[i]
    
