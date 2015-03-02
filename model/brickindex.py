from astropy import wcs
import numpy

class BrickIndex(object):
    def __init__(self, hdudata):
        """ 
            indexing bricks from bricks.fits 
            bricks = fits.open('bricks.fits')
            bi = BrickIndex(bricks[1].data) 
        """

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
        header = dict(
            CTYPE1  = 'RA---TAN',#           / TANgent plane                                  
            CTYPE2  = 'DEC--TAN', #           / TANgent plane                                  
            CRPIX1  =               1800.5, # / Reference x                                    
            CRPIX2  =               1800.5, # / Reference y                                    
            CD1_1   = -7.27777777777778E-05, # / CD matrix                                     
            CD1_2   =                   0., # / CD matrix                                      
            CD2_1   =                   0., # / CD matrix                                      
            CD2_2   = 7.27777777777778E-05, # / CD matrix    
        )
        ra = self.hdudata['RA'][index]  # watch out
        dec = self.hdudata['DEC'][index] # watch out
        header['CRVAL1']  =     ra, # / Reference RA                                   
        header['CRVAL2']  =     dec, # / Reference Dec                                  
        q = wcs.WCS(header)
        return Brick(self.hdudata['BRICKID'][index], 
                self.hdudata['BRICKNAME'][index], wcs)

    def query(self, coord):
        """ 
            querying the indices for given RA and DEC
            RA, DEC = coord

            coord shall be in the same coordinate system of the bricks!
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
        return bid

    def optimize(self, coord):
        """ optimize the ordering of ra, dec,

            This will sort ra, dec by their brickid.

            return sorted_ra, sorted_dec, invert_arg

            invariance:
                sorted_ra[inverg_arg] == ra
                sorted_dec[inverg_arg] == dec
        """
        ra, dec = coord
        bid = self.query_brick(ra, dec)
        arg = bid.argsort()

        invarg = numpy.empty_like(arg)
        invarg[arg] = numpy.arange(len(arg), dtype='i8')
        return ra[arg], dec[arg], invarg


