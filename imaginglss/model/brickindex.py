"""
Python code to provide a high level table/list/index of
bricks, and methods for organizing and querying them and
returning brick objects.

"""

import numpy

from .brick import Brick

__author__ = "Yu Feng and Martin White"
__version__ = "1.0"
__email__  = "yfeng1@berkeley.edu or mjwhite@lbl.gov"



class BrickIndex(object):
    """
    BrickIndex manages a high level table/list/index of bricks.

    It contains methods for organizing and querying bricks
    or returning brick objects.

    The FITS file 'bricks.fits' within each imaging data release
    contains meta-data about each of the bricks in the release.

    We generate an index from bricks.fits , e.g.

    >>> brickdata = fits.read_table('bricks.fits')
    >>> bi = BrickIndex(brickdata) 

    However, a BrickIndex object is usually automatically handled by 
    :py:class:`~model.datarelease.DataRelease`, and one shall use

    >>> dr = DataRelease()
    >>> dr.brickindex

    Notes
    -----
    The bricks in DECALS are first divided as uniform DEC rows. 
    The number of bricks changes with DEC, such that each brick
    as apporximately same area.

    """
    def __init__(self, brickdata):
        """ 
        Parameters
        ----------
        brickdata : array_like
            contents of `bricks.fits`.
    
        """
        self.brickdata = numpy.array(brickdata[:], copy=True)
        self.init_from_state()    

    def __getstate__(self):
        return dict(brickdata=self.brickdata)

    def __setstate__(self, state):

        self.__dict__.update(state)
        self.init_from_state()

    def init_from_state(self):
        brickdata = self.brickdata
        self.ncols = numpy.bincount(brickdata['BRICKROW'])

        self.ROWMAX = brickdata['BRICKROW'].max()
        self.COLMAX = brickdata['BRICKCOL'].max() 

        # fast querying from row col
        self.hash = brickdata['BRICKROW'] * (self.COLMAX + 1) +\
                    brickdata['BRICKCOL']

        self.RA = self.brickdata['RA']
        self.DEC = self.brickdata['DEC']
        self.RA1 = self.brickdata['RA1']
        self.DEC1 = self.brickdata['DEC1']
        self.RA2 = self.brickdata['RA2']
        self.DEC2 = self.brickdata['DEC2']

        self.cache = {}

        self.names_sortarg = self.brickdata['BRICKNAME'].argsort()
        self.names_sorted = self.brickdata['BRICKNAME'][self.names_sortarg]

        assert (self.brickdata['BRICKID'] == numpy.arange(len(self.brickdata)) + 1).all()

    def build(self, NROWS, NCOLS):
        height = 180. / NROWS

        gDEC = numpy.linspace(-90, 90, NROWS + 1, endpoint=True)
        nCOL = numpy.cos(gDEC / 180. * numpy.pi) * NCOLS

        for i, d in enumerate(gDEC):
            pass            
        RA, DEC = coord
        RA = numpy.asarray(RA)
        DEC = numpy.asarray(DEC)
        row = numpy.int32(numpy.floor(DEC * self.ROWMAX / 180 + 360. + 0.5))
        row = numpy.clip(row, 0, self.ROWMAX)
        ncols = self.ncols[row]
        col = numpy.int32(numpy.floor(RA * ncols / 360. ))
        hash = row * (self.COLMAX + 1) + col
        ind = self.hash.searchsorted(hash)
         
    def __len__(self):
        return len(self.brickdata)

    def get_brick(self, index):
        """
        Obtain a single brick from index.

        The returned Brick object is immutable; there is only
        a single instance per BrickIndex.
    
        Parameters
        ----------
        index : integer
            The internal index of a brick. This differ from BRICKID
            in the catalogue.

        Returns
        -------
        brick : :py:class:`~model.brick.Brick`
            A brick object at index.

        Notes
        -----
        The lack of a vector version is on purpose to emphasize
        we are dealing with objects. 

        Use :py:meth:`get_bricks` to create a list of objects.

        """
        if index not in self.cache:
            brick = Brick(
                    index,
                    self.brickdata['BRICKID'][index], 
                    self.brickdata['BRICKNAME'][index], 
                    self.brickdata['RA'][index], 
                    self.brickdata['DEC'][index], 
                    self.brickdata['RA1'][index], 
                    self.brickdata['RA2'][index], 
                    self.brickdata['DEC1'][index], 
                    self.brickdata['DEC2'][index])
            self.cache[index] = brick
        return self.cache[index]

    def get_bricks(self, indices):
        """
        Obtain a list of bricks.

        Parameters
        ----------
        indices : array_like
            Indices to obtain bricks

        Returns
        -------
        bricks  : list
            A list of bricks.

        """ 
        return [self.get_brick(i) for i in indices]

    def search_by_name(self, brickname):
        """
        Search the brickindex for bricks with a given name.

        Parameters
        ----------
        brickname : string
            the matching BRICKNAME.

        Returns
        -------
        index  : integer
            The internal index of bricks.

        Notes
        -----
        This function is not vectorized.
        """
        foo = numpy.empty(1, self.brickdata['BRICKNAME'].dtype)
        foo[0] = brickname
        bid = self.names_sortarg[self.names_sorted.searchsorted(foo[0])]
        return bid
 
    def search_by_id(self, brickid):
        """
        Search the brickindex for bricks with a given brickid.
        
        Parameters
        ----------
        brickid : integer or array_like
            the matching BRICKID or list of BRICKIDs

        Returns
        ------- 
        index   : integer or array_like
            the internal index of brick or bricks (for array_like input)

        Notes
        -----
        Get Brick objects with :py:meth:`get_bricks` or :py:meth:`get_brick`.
        
        """
        rt = self.brickdata['BRICKID'].searchsorted(brickid)
        rt = rt.clip(0, len(self.brickdata) - 1)
        found = self.brickdata['BRICKID'][rt] == brickid
        if not found.all():
            raise IndexError("Some brickid %s are not found"%\
              str(brickid[~found]))
        return rt

    def query_internal(self, coord):
        """ 
        Returns the internal index of a brick at coord=(RA,DEC) in
        decimal degrees.  

        Parameters
        ----------
        coord  : array_like
            coord = (RA, DEC) in degrees, vectorized.

        Notes
        -----
        Get Brick objects with :py:meth:`get_bricks` or :py:meth:`get_brick`.

        """
        RA, DEC = coord
        RA = numpy.asarray(RA)
        DEC = numpy.asarray(DEC)
        row = numpy.int32(numpy.floor(DEC * self.ROWMAX / 180 + 360. + 0.5))
        row = numpy.clip(row, 0, self.ROWMAX)
        ncols = self.ncols[row]
        col = numpy.int32(numpy.floor(RA * ncols / 360. ))
        hash = row * (self.COLMAX + 1) + col
        ind = self.hash.searchsorted(hash)
        return ind

    def optimize(self, coord, return_index=False, return_inverse=False):
        """
        Optimize the ordering of coord=(RA,DEC) to make later queries faster.

        RA and DEC are sorted by their brickid, to group future queries.

        Parameters
        ----------
        return_inverse : boolean
            if True, returns the array that can be used to
            reconstruct coord from the returned coord.
        return_index   : boolean
            if True, returns the array that can be used to
            construct sorted_coord from coord:

        Returns
        -------
        Depending or return_inverse and return_index, may return
        (sorted_coord, indices, iindices), (sorted_coord, iindices)
        (sorted_coord, indices), or sorted_coord

        sorted_coord : array_like
            Optimized coord array, that is sorted by bricks.
        iindeces     : array_like
            sorted_ra[ iindices] == ra
            sorted_dec[iindices] == dec
        indices      : array_like
            sorted_ra == ra[indices]
            sorted_dec == dec[indices]

        """
        coord = numpy.array(coord)
        bid = self.query_internal(coord)
        arg = bid.argsort()
        #
        if return_inverse:
            invarg = numpy.empty_like(arg)
            invarg[arg] = numpy.arange(len(arg), dtype='i8')
            if return_index:
                return numpy.array(coord[:, arg]), arg, invarg
            else:
                return numpy.array(coord[:, arg]), invarg
        else:
            if return_index:
                return numpy.array(coord[:, arg]), arg
            else:
                return numpy.array(coord[:, arg])


