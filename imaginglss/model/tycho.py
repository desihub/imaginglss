import numpy
from ..utils import fits
from kdcount import KDTree

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

class Tycho(numpy.ndarray):
    """ Representing a Tycho catalogue that
        is used to veto objects near stars 
    """
    def __new__(kls, path):
        self = fits.read_table(path).view(type=kls)
        return self

    def __init__(self, path):
        self.bmag = self['BMAG']
        self.vmag = self['VMAG']
        self.varflag = self['VARFLAG']
        pos = radec2pos(self['RA'], self['DEC'])
        self.tree = KDTree(pos)

    def veto(self, coord, vetotype):
        """
            Returns a veto mask for coord.

            Parameters
            ----------
            coord : (RA, DEC)
            vetotype : a function that f(bmag, vmag) -> Radius in degrees

            Returns
            -------
            Vetomask : True for veto, False for keep.

        """
        R = vetotype(self.bmag, self.vmag) # in degrees

        # convert to euclidean distance
        R = 2 * numpy.sin(numpy.radians(R) * 0.5)

        Rmax = R.max()

        pos = radec2pos(coord[0], coord[1])
        other = KDTree(pos)
        vetoflag = numpy.zeros(len(pos), dtype='?')

        def process(r, i, j):
            # i is tycho, j is objects
            rcut = R[i]
            jcut = j[r < rcut]
            vetoflag[jcut] |= True

        self.tree.root.enum(other.root, Rmax, process)
        return vetoflag
