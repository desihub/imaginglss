import numpy

from astropy.io import fits
from brick import Brick

class BrickIndex(object):
    def __init__(self, hdudata):
        """ 
            indexing bricks from bricks.fits 
            bricks = fits.open('bricks.fits')
            bi = BrickIndex(bricks[1].data) 

            hdudata can also be a file name.
        """
        if isinstance(hdudata, basestring):
            hdudata = fits.open(hdudata)[1].data[:]
            
        self.hdudata = numpy.array(hdudata[:], copy=True)
        self.ncols = numpy.bincount(hdudata['BRICKROW'])

        self.ROWMAX = hdudata['BRICKROW'].max()
        self.COLMAX = hdudata['BRICKCOL'].max() 

        # fast querying from row col
        self.hash = hdudata['BRICKROW'] * (self.COLMAX + 1) + hdudata['BRICKCOL']

        assert (self.hdudata['BRICKID'] == numpy.arange(len(self.hdudata)) + 1).all()

    def __getitem__(self, index):
        """ create a single brick from bid """
        # template header to feed wcs
        # fill CRVAL1, CRVAL2 later
        return Brick(self.hdudata['BRICKID'][index], 
                self.hdudata['BRICKNAME'][index], 
                self.hdudata['RA'][index], 
                self.hdudata['DEC'][index], 
                self.hdudata['RA1'][index], 
                self.hdudata['RA2'][index], 
                self.hdudata['DEC1'][index], 
                self.hdudata['DEC2'][index])

    def query(self, coord):
        """ 
            querying the indices for given RA and DEC
            RA, DEC = coord

            coord shall be in the same coordinate system of the bricks!
            returns the internal index of the bricks.
            get real Brick objects by iterating over the result and use brickindex[i]
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

    def optimize(self, coord):
        """ optimize the ordering of ra, dec,

            This will sort ra, dec by their brickid.

            return sorted_ra, sorted_dec, invert_arg

            invariance:
                sorted_ra[inverg_arg] == ra
                sorted_dec[inverg_arg] == dec
        """
        coord = numpy.array(coord)
        bid = self.query(coord)
        arg = bid.argsort()

        invarg = numpy.empty_like(arg)
        invarg[arg] = numpy.arange(len(arg), dtype='i8')
        return numpy.array(coord[:, arg]), invarg


