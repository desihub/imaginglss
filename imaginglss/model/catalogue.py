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
        if isinstance(catalogue, CachedCatalogue):
            return catalogue[self.name][:]
        else:
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

    COLUMNS = [
        'BRICK_PRIMARY',
        'RA',
        'DEC',
        'DECAM_FLUX_IVAR',
        'DECAM_MW_TRANSMISSION',
        'DECAM_PSFSIZE',
        'DECAM_NOBS',
        'DECAM_ANYMASK',
        'DECAM_DEPTH',
        'DECAM_FLUX',
        'WISE_FLUX',
        'WISE_FLUX_IVAR',
        'WISE_MW_TRANSMISSION',
        'TYPE',
        'SHAPEDEV_R',
        'SHAPEEXP_R',
    ]

    def __init__(self, bricks, format_filename, aliases):

        filenames = [ format_filename(brick) for brick in bricks]
        bricknames = [ brick.name for brick in bricks]

        self.filenames = dict(zip(bricknames, filenames))

        self.aliases = dict([(new, (old, transform)) 
                for old, new, transform in aliases])

        data = []
        for brick in bricks:
            data.append(self.open(brick))

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


class CachedCatalogue(ColumnStore):
    """
    Class for handling object catalogs associated with a data release.

    Access via the :py:attr:`catalogue` attribute of :py:class:`DataRelease`
    object.

    Notes
    -----
    The catalogs are contained in many small FITS files. Accesing them 
    directly is slow. 
    The columns must be first converted from the many-small file original
    format to a cache format via a script, :code:`scripts/build_cache.py`.
    The scripts calls the function :py:meth:`build_cache` to build the cache
    for chunks of catalogue in parallel, then write them to the correct cache
    directory.

    This class caches the information on disk for speed. 
    Only columns that are accessed are loaded into memory.

    Examples
    --------
    >>> d = DECALS()
    >>> cat = d.datarelease.catalogue
    >>> print cat['BRICK_PRIMARY'][:100]
    >>> print cat['BRICK_PRIMARY'][10]
    >>> print cat['BRICK_PRIMARY'][:] # may be huge!

    Parameters
    ----------
    cachedir  : string
        the location for caching.
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
    filenames : dict
        A mapping between the brickname and the filenames.
    
    """

    def __init__(self, cachedir, bricks, format_filename, aliases):

        filenames = [ format_filename(brick) for brick in bricks]
        bricknames = [ brick.name for brick in bricks]

        self.filenames = dict(zip(bricknames, filenames))

        self.aliases = dict([(new, (old, transform)) 
                for old, new, transform in aliases])
        self.cachedir = cachedir

        try:
            self._size = filehandler.size(self.cachedir, 'BRICK_PRIMARY')
            columns = filehandler.list(self.cachedir)
            self._dtype = numpy.dtype([(key, columns[key]) for key in columns])
        except:
            self._size = None
            self._dtype = None
        ColumnStore.__init__(self)

    @property
    def size(self):
        return self._size

    @property
    def dtype(self):
        return self._dtype

    def open(self, brick):
        return fits.read_table(self.filenames[brick.name])

    def build_cache(self, filenames):
        """
        Build cache of the catalogue files listed in filenames

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
       
        Returns
        -------
        concatenated catalogue for filenames

        """ 
        
        cachedir = self.cachedir

        # get the first filename
        dtype = None
        for i, fn in enumerate(self.filenames.values()):
            try:
                first = fits.read_table(fn)
                dtype = subdtype(uppercase_dtype(first.dtype), Catalogue.COLUMNS)
            except:
                if i < 10 and i != len(self.filenames) - 1:
                    pass
                else :
                    raise
                
        total = len(filenames)

        data = numpy.empty(0, dtype)

        for filename in filenames:
            table = fits.read_table(filename)
            table = table.view(uppercase_dtype(table.dtype))
            data1 = numpy.zeros(len(table), dtype)
            for name in table.dtype.names:
                # only preserve those in both 'first' and all
                if name not in dtype.names: continue
                data1[name][...] = table[name]

            data = numpy.append(data, data1)

        return data

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
        if isinstance(column, basestring) and column in self.aliases:
            old, transform = self.aliases[column]
            return TransformedColumn(self, old, transform)
        else:
            return ColumnStore.__getitem__(self, column)

    def fetch(self, column, start, end):
        if not self.check_cache():
            raise CacheExpired("The cache is too old. Regenerate it with scripts/build_cache.py")

        return filehandler.read(self.cachedir, [column], offset=start, count=end-start)[column]

    def __repr__(self):
        return 'CachedCatalogue: %s' % str(self.dtype)

    def neighbours(self, coord, sep):
        pass

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
