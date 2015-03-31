# Python class for handling object catalogs associated with
# a data release.  The catalogs are obtained from FITS files.
# This class does some caching for speed.

import numpy

from utils import fits
from utils.columnstore import DiskColumnStore
import multiprocessing.pool

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

def uppercase_dtype(dtype):
    pairs = dict([(key.upper(), dtype.fields[key]) for key in dtype.names])
    dtype = numpy.dtype(pairs)
    return dtype

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
    def __init__(self, cachedir, filenames, aliases):
        """
        cachedir is the location for caching.
        filenames are a list of fits file names for the DR
        aliases is a list of fields to rename (fields in old DRs are renamed to newer DRs,
        so that code written for new DR can work on old DRs.
        """
        self.filenames = filenames
        fn = filenames[0]
        first = fits.read_table(fn)
        self.aliases = dict([(new, (old, transform)) 
                for old, new, transform in aliases])
        dtype = uppercase_dtype(first.dtype)
        DiskColumnStore.__init__(self, cachedir, dtype)

    def __getitem__(self, column):
        if column in self.aliases:
            old, transform = self.aliases[column]
            return transform(self[old])
        else:
            return DiskColumnStore.__getitem__(self, column)

    def fetch(self, column):
        def readafile(fn):
            # FIXME: this reads in the full table .. maybe not a good idea
            # but the files are small ..
            data = numpy.array(fits.read_table(fn))
            data = data.view(dtype=uppercase_dtype(data.dtype))
            return numpy.array(data[column], dtype=self.dtype[column].base)
        pool = multiprocessing.pool.ThreadPool()
            #cat = [ readafile(fn) for fn in self.filenames]
        cat = pool.map(readafile, self.filenames)
        cat = numpy.concatenate(cat, axis=0)
        return cat

    def __repr__(self):
        return 'Catalogue: %s' % str(self.dtype)

    def neighbours(self, coord, sep):
        pass
