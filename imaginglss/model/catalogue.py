"""
Python class for handling object catalogs associated with
a data release.  The catalogs are obtained from FITS files.
This class does some caching for speed.

"""
import numpy

from ..utils import fits
from ..utils import filehandler
from ..utils.columnstore import ColumnStore
from ..utils.npyquery import Column as ColumnBase

try:
  basestring
except NameError:
  basestring = str

__author__ = "Yu Feng and Martin White"
__version__ = "1.0"
__email__  = "yfeng1@berkeley.edu or mjwhite@lbl.gov"

class C(ColumnBase):
    """ `C` provides a shorthand for querying the catalogue
        with the :py:mod:`imaginglss.utils.npyquery` mini language 
     """
    def visit(self, catalogue):
        # Rows
        return catalogue[self.name]

def coord2xyz(coord):
    """
    Given coord=(RA,DEC) returns unit vectors, nhat.  A helper function.
    
    Parameters
    ----------
    coord  : array_like
        coord = (RA, DEC) in degrees.
    
    Returns
    -------
    vector : array_like
        Unit vectors corresponding to RA, DEC, in (, 3).

    """
    RA, DEC = coord
    xyz = numpy.empty(len(RA), ('f4', 3))
    c = numpy.cos(DEC / 180. * numpy.pi)
    xyz[:, 0] = c * numpy.sin(RA / 180. * numpy.pi)
    xyz[:, 1] = c * numpy.cos(RA / 180. * numpy.pi)
    xyz[:, 2] = numpy.sin(DEC / 180. * numpy.pi)
    return xyz.T

def uppercase_dtype(dtype):
    """ Convert a dtype to upper case. A helper function.
        
        Do not use.
    """
    pairs = dict([(key.upper(), dtype.fields[key]) for key in dtype.names])
    dtype = numpy.dtype(pairs)
    return dtype

def subdtype(dtype, columns):
    try:
        return numpy.dtype([(key, dtype[key]) for key in columns])
    except KeyError as e:
        raise KeyError("%s : candidates are %s." %
                (str(e), str(sorted(dtype.names))))
        
def native_dtype(dtype):
    """ Convert a dtype to native dtype. A helper function.
        
        Do not use.
    """
    return dtype.newbyteorder('=')

class CacheExpired(RuntimeError):
    pass

class TransformedColumn(object):
    def __init__(self, ref, columns, transform):
        if not isinstance(columns, (tuple, list)):
            columns = [columns]
        self.ref = ref
        self.columns = columns
        self.transform = transform

    def __getitem__(self, index):
        args = tuple([  self.ref[c][index]
                        for c in self.columns])
        return self.transform(*args)

class Catalogue(object):
    """
    Parameters
    ----------
    bricks: list
        a list of bricks names that the catalogue covers.
    format_filename : function
        a function converts a brick object to a filename of the tractor
        catalogue
    aliases   : list
        a list of fields to transform; this is to support migration
        of schema from older data release to newer ones. The list
        is of from (oldname, newname, transformfunction)

    Attributes
    ----------
    dtype : dtype
        A container of the data type of columns
        in :py:class:`numpy.dtype`
    """

    def __init__(self, bricks, format_filename, aliases, columns):

        filenames = [ format_filename(brick) for brick in bricks]
        bricknames = [ brick.name for brick in bricks]

        self.filenames = dict(zip(bricknames, filenames))

        self.aliases = dict([(new, (old, transform)) 
                for old, new, transform in aliases])

        data = []
        for brick in bricks:
            data.append(self.open(brick))

        self.COLUMNS = columns
        self.data = numpy.concatenate(data)

    @property
    def size(self):
        return len(self.data)

    def __len__(self):
        return len(self.data)

    @property
    def dtype(self):
        return self.data.dtype

    def open(self, brick):
        data = fits.read_table(self.filenames[brick.name])
        dtype = [(column, data.dtype[column]) for column in self.COLUMNS]
        data_compressed = numpy.empty(shape=data.shape, dtype=dtype)
        for column in self.COLUMNS:
            data_compressed[column][...] = data[column]
        return data_compressed

    def __getitem__(self, column):
        if isinstance(column, basestring) and column in self.aliases:
            old, transform = self.aliases[column]
            return TransformedColumn(self, old, transform)
        else:
            return self.data[column]

    def __repr__(self):
        return 'Catalogue: %s' % str(self.dtype)

class BigFileCatalogue(ColumnStore):
    """
    """

    def __init__(self, cachedir, aliases):
        import bigfile
        self.cachedir = cachedir

        with bigfile.BigFile(cachedir, create=True) as bf:
            bd = bigfile.BigData(bf)

            self._size = bd.size
            self._dtype = bd.dtype

        self.aliases = dict([(new, (old, transform)) 
                for old, new, transform in aliases])
        ColumnStore.__init__(self)

    @property
    def size(self):
        return self._size

    @property
    def dtype(self):
        return self._dtype

    def open(self, brick):
        raise RuntimeError("FIXME: currently cannot open a brick from a sweep.")

    def __getitem__(self, column):
        print(self.aliases)
        if isinstance(column, basestring) and column in self.aliases:
            old, transform = self.aliases[column]
            return TransformedColumn(self, old, transform)
        else:
            return ColumnStore.__getitem__(self, column)

    def fetch(self, column, start, end):
        import bigfile
        with bigfile.BigFile(self.cachedir) as bf:
            return bf[column][start:end]

    def __repr__(self):
        return 'BigFileCatalogue: %s' % str(self.dtype)

    def neighbours(self, coord, sep):
        pass
