from astropy import wcs
from astropy.io import fits
import numpy

class BrickIndex(object):
    def __init__(self, hdudata):
        """ 
            indexing bricks from bricks.fits 
            bricks = fits.open('bricks.fits')
            bi = BrickIndex(bricks[1].data) 
            
        """

        #    FIXME:
        #    hard coded numbers: 10000 (max number of rows per col)
        self.hdudata = hdudata[:].copy()
        self.ncols = numpy.bincount(hdudata['BRICKROW'])
        self.hash = hdudata['BRICKROW'] * 10000 + hdudata['BRICKCOL']
        assert (self.hdudata['BRICKID'] == numpy.arange(len(self.hdudata)) + 1).all()

    def query_brick(self, RA, DEC):
        """ 
            finds the brick index (BRICKID - 1) for given RA and DEC 

        """
        #    FIXME:
        #    hard coded number: 720(~number of cols - 1), 10000 (max number of rows per col)
        #    4, 360. 0.5 the spacing of cols (found via a poly fit, clear this up!)
        RA = numpy.asarray(RA)
        DEC = numpy.asarray(DEC)
        row = numpy.int32(numpy.floor(DEC * 4 + 360. + 0.5))
        row = numpy.clip(row, 0, 720)
        ncols = self.ncols[row]
        col = numpy.int32(numpy.floor(RA * ncols / 360. ))
        hash = row * 10000 + col
        ind = self.hash.searchsorted(hash)
        return ind

    def query(self, RA, DEC):
        """ 
            This will return the brickid and pix x, y (0, 0 as origin)

            RA DEC must be arrays of same length (no scalars yet!)

            We can make it faster by packing calls to wcs.all_world2pix of
            the same brick together. 
        """
        brk = self.query_brick(RA, DEC)
        
#        ubrk, ind = numpy.unique(brk, return_inverse=True)
        pix = numpy.empty((len(brk), 3), dtype='i4') 
        pix[:, 0] = self.hdudata['BRICKID'][brk]

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

        # now loop and query!
        oldbrk = None
        for i in range(len(brk)):
            if brk[i] != oldbrk:
                ra = self.hdudata['RA'][brk[i]]
                dec = self.hdudata['DEC'][brk[i]]
                header['CRVAL1']  =     ra, # / Reference RA                                   
                header['CRVAL2']  =     dec, # / Reference Dec                                  

                q = wcs.WCS(header)
                oldbrk = brk[i]
             
            dat = numpy.array([[RA[i], DEC[i]]])
            r = q.all_world2pix(dat, 0)[0]
            pix[i, 1] = r[0]
            pix[i, 2] = r[1]

        return pix

    def test(self):
        """ no testing on RA; this makes sure the DEC and centers edges are correctly handled """
        rows = self.hdudata[self.query_brick(self.hdudata['RA'], self.hdudata['DEC'])]
        print 'failed', (rows['BRICKID'] != self.hdudata['BRICKID']).nonzero()
        rows = self.hdudata[self.query_brick(self.hdudata['RA'], self.hdudata['DEC'] - 0.125)]
        print 'failed', (rows['BRICKID'] != self.hdudata['BRICKID']).nonzero()
        rows = self.hdudata[self.query_brick(self.hdudata['RA'], self.hdudata['DEC'] + 0.124)]
        print 'failed', (rows['BRICKID'] != self.hdudata['BRICKID']).nonzero()

def load(repo, brickid, x, y):
    """ repo is a string with %(brickid),

        load all pixels indexed by brickid, x, y from the primary HDU of fits files in repo,

        eg:
            bxy = bi.query([243.6] * 6, numpy.arange(11.75 - 0.12, 11.75 + 0.12, 0.04))
            print load('coadd/depth-%(brickid)d-z.fits.gz', *bxy.T)

        currently we avoid reopening the files if the brickid is continous.
        This can be done better!
    """
    pixels = numpy.empty(len(brickid))
    oldid = None
    for i in range(len(brickid)):
        if brickid[i] != oldid:
            image = fits.open(repo % dict(brickid=brickid[i]))[0].data
            oldid = brickid[i]
        pixels[i] =  image[(x[i], y[i])]
    return pixels

if __name__ == '__main__':
    bricks = fits.open('bricks.fits')
    bi = BrickIndex(bricks[1].data) 
    print bricks[1].data[398599 - 1]
    #print bricks[1].data.dtype
    #print bricks[1].data[900]
    bxy = bi.query([243.6] * 6, numpy.arange(11.75 - 0.12, 11.75 + 0.12, 0.04))
    print bxy
    print load('coadd/depth-%(brickid)d-z.fits.gz', *bxy.T)
    #bi.test()

