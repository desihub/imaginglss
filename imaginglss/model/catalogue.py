"""
Python class for handling object catalogs associated with
a data release.  The catalogs are obtained from FITS files.
This class does some caching for speed.

"""
import numpy

from ..utils import fits
from ..utils import filehandler
from ..utils.columnstore import ColumnStore

__author__ = "Yu Feng and Martin White"
__version__ = "1.0"
__email__  = "yfeng1@berkeley.edu or mjwhite@lbl.gov"




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

def native_dtype(dtype):
    """ Convert a dtype to native dtype. A helper function.
        
        Do not use.
    """
    return dtype.newbyteorder('=')

class CacheExpired(RuntimeError):
    pass

class Catalogue(ColumnStore):
    """
    Class for handling object catalogs associated with a data release.

    The catalogs are contained in many small FITS files. Accesing them 
    directly is slow. 
    The columns must be first converted from the many-small file original
    format to a cache format via :py:func:`build_cache`.

    This class caches the information on disk for speed. 
    Only columns that are accessed are loaded into memory.

    Parameters
    ----------
    cachedir  : string
        the location for caching.
    filenames : list
        a list of fits file names that the catalogue is stored.
    aliases   : list
        a list of fields to transform; this is to support migration
        of schema from older data release to newer ones. The list
        is of from (oldname, newname, transformfunction)

    """
    def __init__(self, cachedir, filenames, aliases):
        self.filenames = filenames
        fn = filenames[0]
        first = fits.read_table(fn)
        self.aliases = dict([(new, (old, transform)) 
                for old, new, transform in aliases])
        dtype = native_dtype(uppercase_dtype(first.dtype))

        self.cachedir = cachedir
        ColumnStore.__init__(self, dtype)

    def build_cache(self):
        """
        Build Cache of the catalogue.

        The fits files are converted to file handler format.
        Each column becomes a single file.

        Notes
        -----
        This shall be run before using the catalogue. 
        And it takes a long time.

        Parameters
        ----------
        cachedir : string
            directory for holding the cache
        filenames : list
            list of FITS files to read from
       
        """ 
        
        filenames = self.filenames
        cachedir = self.cachedir

        fn = filenames[0]
        first = fits.read_table(fn)
        dtype = uppercase_dtype(first.dtype)

        filehandler.write(cachedir, first.view(dtype), mode='w')

        for filename in filenames[1:]:
            table = fits.read_table(filename)
            filehandler.write(cachedir, table.view(dtype), mode='a')

        d = {
            'nfiles': numpy.array([len(filenames)],dtype='i8') 
            }

        filehandler.write(cachedir, d)

    def check_cache(self):
        """ Check if cache is consistent 
            
            Returns
            -------
            consistent : boolean
                True if consistent, Falst if need to rebuild with :py:meth:`build_cache`.
        """
        Nfile = filehandler.read(self.cachedir, ["nfiles"])['nfiles'][0]
        try:
            Nfile = filehandler.read(self.cachedir, ["nfiles"])['nfiles'][0]
        except filehandler.MissingColumn:
            Nfile = -1
        return Nfile == len(self.filenames)
        
    def __getitem__(self, column):
        if column in self.aliases:
            old, transform = self.aliases[column]
            return transform(self[old])
        else:
            return ColumnStore.__getitem__(self, column)

    def fetch(self, column):
        if not self.check_cache():
            raise CacheExpired("The cache is too old. Regenerate it with imaginglss.model.catalogue.build_cache")

        return filehandler.read(self.cachedir, [column])[column]

    def __repr__(self):
        return 'Catalogue: %s' % str(self.dtype)

    def neighbours(self, coord, sep):
        pass
