# Python class for handling object catalogs associated with
# a data release.  The catalogs are obtained from FITS files.
# This class does some caching for speed.

import numpy

from utils import fits
from columnstore import DiskColumnStore


__author__ = "Yu Feng and Martin White"
__version__ = "1.0"
__email__  = "yfeng1@berkeley.edu or mjwhite@lbl.gov"




def coord2xyz(coord):
    """
    Given coord=(RA,DEC) returns unit vectors, nhat.  A helper function.
    """
    RA, DEC = coord
    xyz = numpy.empty(len(RA), ('f4', 3))
    c = numpy.cos(DEC / 180. * numpy.pi)
    xyz[:, 0] = c * numpy.sin(RA / 180. * numpy.pi)
    xyz[:, 1] = c * numpy.cos(RA / 180. * numpy.pi)
    xyz[:, 2] = numpy.sin(DEC / 180. * numpy.pi)
    return xyz.T

class Catalogue(DiskColumnStore):
    """
    Class for handling object catalogs associated with a data release.
    The catalogs are contained in FITS files, but this class caches the
    information for speed and only columns that are accessed are loaded
    into memory.
    The columns are cached into cachedir, which is initialized 
    by the DataRelease object.
    We assume the catalogue is split into multiple files.
    and the list of file names are provided as an input.
    The main method is fetch.
    """
    def __init__(self, cachedir, filenames):
        """
        cachedir is the location for caching.
        filenames are a list of fits file names for the DR
        """
        self.filenames = filenames
        fn = filenames[0]
        first = fits.read_table(fn)
        DiskColumnStore.__init__(self, cachedir, first.dtype)

    def fetch(self, column):
        cat = [numpy.array(fits.read_table(fn), copy=True)[column]
            for fn in self.filenames]
        cat = numpy.concatenate(cat, axis=0)
        return cat

    def __repr__(self):
        return 'Catalogue: %s' % str(self.dtype)

    def neighbours(self, coord, sep):
        pass
