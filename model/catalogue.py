import numpy
from astropy.io import fits

def coord2xyz(coord):
    RA, DEC = coord
    xyz = numpy.empty(len(RA), ('f4', 3))
    c = numpy.cos(DEC / 180. * numpy.pi)
    xyz[:, 0] = c * numpy.sin(RA / 180. * numpy.pi)
    xyz[:, 1] = c * numpy.cos(RA / 180. * numpy.pi)
    xyz[:, 2] = numpy.sin(DEC / 180. * numpy.pi)
    return xyz.T

class Catalogue(object):
    """ Catalogue. 
        Only columns that are accessed are loaded to the memory.

        This shall be split into a ColumnStore interface, with optional
        on-disk caching facilities.
    """
    def __init__(self, filenames):
        self.filenames = filenames
        self.cache = {}
        fn = filenames[0]
        first = fits.open(fn)[1].data
        self.dtype = first.dtype

    def __getitem__(self, column):
        if column not in self.cache:
            cat = [numpy.array(fits.open(fn)[1].data[column], copy=True)
                for fn in self.filenames]
            cat = numpy.concatenate(cat)
            self.cache[column] = cat
        return self.cache[column] 

    def __contains__(self, column):
        return column in self.dtype.names

    def __setitem__(self, column):
        raise "Unsupported"

    def forget(self, column):
        del self.cache[column]

    def __repr__(self):
        return 'Catalogue: %s' % str(self.dtype)

    def neighbours(self, coord, sep):
        pass
