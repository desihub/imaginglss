# Python code to provide a high level table/list/index of
# bricks, and methods for organizing and querying them and
# returning brick objects.
#
#

import numpy

from brick import Brick

__author__ = "Yu Feng and Martin White"
__version__ = "1.0"
__email__  = "yfeng1@berkeley.edu or mjwhite@lbl.gov"



class BrickIndex(object):
    """
    Manages a high level table/list/index of bricks.
    Contains methods for organizing and querying bricks
    or returning brick objects.
    """
    def __init__(self, brickdata):
        """ 
        The FITS file 'bricks.fits' within each imaging data release
        contains meta-data about each of the bricks in the release.

        We generate an index from bricks.fits , e.g.
        brickdata = fits.read_table('bricks.fits')
        bi = BrickIndex(brickdata) 
    
        """
        self.brickdata = numpy.array(brickdata[:], copy=True)

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
         
    def __len__(self):
        return len(self.brickdata)

    def get_brick(self, index):
        """
        Obtain a single brick from index.
        The returned Brick object is immutable; there is only
        a single instance per BrickIndex.
        The lack of a vector version is on purpose to emphasize
        we are dealing with objects. (YF: maybe not a good idea?)
        """
        if index not in self.cache:
            brick = Brick(self.brickdata['BRICKID'][index], 
                    self.brickdata['BRICKNAME'][index], 
                    self.brickdata['RA'][index], 
                    self.brickdata['DEC'][index], 
                    self.brickdata['RA1'][index], 
                    self.brickdata['RA2'][index], 
                    self.brickdata['DEC1'][index], 
                    self.brickdata['DEC2'][index])
            self.cache[index] = brick
        return self.cache[index]

    def search_by_name(self, brickname):
        """
        search the brickindex for bricks with a given name.
        returns the internal index of bricks.
        Get Brick objects by iterating over the result and using get_brick(i)
        """
        foo = numpy.empty(1, self.brickdata['BRICKNAME'].dtype)
        foo[0] = brickname
        bid = self.names_sortarg[self.names_sorted.searchsorted(foo[0])]
        return bid
 
    def search_by_id(self, brickid):
        """
        search the brickindex for bricks with a given brickid.
        returns the internal index of bricks.
        Get Brick objects by iterating over the result and using get_brick(i)
        """
        rt = self.brickdata['BRICKID'].searchsorted(brickid)
        rt = rt.clip(0, len(self.brickdata) - 1)
        found = self.brickdata['BRICKID'][rt] == brickid
        if not found.all():
            raise IndexError("Some brickid %s are not found"%\
              str(brickid[~found]))
        return rt

    def query(self, coord):
        """ 
        Returns the internal index of a brick at coord=(RA,DEC) in
        decimal degrees.  This internal index can differ from BRICKID.
        Get Brick objects by iterating over the result and using get_brick(i)
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

    def optimize(self, coord, return_index=False, return_inverse=True):
        """
        Optimize the ordering of coord=(RA,DEC) to make later queries faster.
        RA and DEC are sorted by their brickid, to group future queries.

        if return_inverse is True, returns the array that can be used to
        reconstruct coord from the returned coord.

                sorted_ra[ iindices] == ra
                sorted_dec[iindices] == dec

        if return_index is True, returns the array that can be used to
        construct sorted_coord from coord:

                sorted_ra == ra[indices]
                sorted_dec == dec[indices]

        returns sorted_coord, indices, iindices or
                sorted_coord, iindices or
                sorted_coord, indices or
                sorted_coord
        """
        coord = numpy.array(coord)
        bid = self.query(coord)
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


