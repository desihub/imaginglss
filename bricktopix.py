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
        self.hdudata = numpy.array(hdudata[:], copy=True)
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

    def revert(self, brickid, x, y):
        """ brickid, x, y -> RA, DEC """
        pix = numpy.empty((len(brickid), 2), dtype='f8') 

        # translate to BRICKID + 1 by searching (could done - 1)
        brickid = self.hdudata['BRICKID'].searchsorted(brickid)

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
        oldbrk = -1 
        start = 0
        q = None
        for i in range(len(brickid) + 1):
            if not (i == len(brickid) or brickid[i] != oldbrk): continue
            sl = slice(start, i)

            if i != 0:
                dat = numpy.array((x[sl], y[sl])).T
                r = q.all_pix2world(dat, 0)
                pix[sl, 0] = r[:, 0]
                pix[sl, 1] = r[:, 1]

            if i != len(brickid):
                #advance
                ra = self.hdudata['RA'][brickid[i]]
                dec = self.hdudata['DEC'][brickid[i]]
                header['CRVAL1']  =     ra, # / Reference RA                                   
                header['CRVAL2']  =     dec, # / Reference Dec                                  

                q = wcs.WCS(header)

                oldbrk = brickid[i]
                start = i

        return pix.T

    def query(self, RA, DEC):
        """ 
            This will return the brickid and pix x, y (0, 0 as origin)

            RA DEC must be arrays of same length (no scalars yet!)

            We can make it faster by packing calls to wcs.all_world2pix of
            the same brick together. 
        """
        brk = self.query_brick(RA, DEC)
        
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
        oldbrk = -1 
        start = 0
        q = None
        for i in range(len(brk) + 1):
            if not (i == len(brk) or brk[i] != oldbrk): continue
            sl = slice(start, i)

            if i != 0:
                dat = numpy.array([RA[sl], DEC[sl]], dtype='f8').T
                r = q.all_world2pix(dat, 0)
                pix[sl, 1] = r[:, 0]
                pix[sl, 2] = r[:, 1]

            if i != len(brk):
                #advance
                ra = self.hdudata['RA'][brk[i]]
                dec = self.hdudata['DEC'][brk[i]]
                header['CRVAL1']  =     ra, # / Reference RA                                   
                header['CRVAL2']  =     dec, # / Reference Dec                                  

                q = wcs.WCS(header)

                oldbrk = brk[i]
                start = i

        return pix.T

    def test(self):
        """ no testing on RA; this makes sure the DEC and centers edges are correctly handled """
        rows = self.hdudata[self.query_brick(self.hdudata['RA'], self.hdudata['DEC'])]
        print 'failed', (rows['BRICKID'] != self.hdudata['BRICKID']).nonzero()
        rows = self.hdudata[self.query_brick(self.hdudata['RA'], self.hdudata['DEC'] - 0.125)]
        print 'failed', (rows['BRICKID'] != self.hdudata['BRICKID']).nonzero()
        rows = self.hdudata[self.query_brick(self.hdudata['RA'], self.hdudata['DEC'] + 0.124)]
        print 'failed', (rows['BRICKID'] != self.hdudata['BRICKID']).nonzero()

def load(repo, brickid, x, y, default=numpy.nan):
    """ repo is a string with %(brickid),

        load all pixels indexed by brickid, x, y from the primary HDU of fits files in repo,

        eg:
            bxy = bi.query([243.6] * 6, numpy.arange(11.75 - 0.12, 11.75 + 0.12, 0.04))
            print load('coadd/depth-%(brickid)d-z.fits.gz', *bxy.T)

        currently we avoid reopening the files if the brickid is continous.
        This can be done better!
    """
    pixels = numpy.empty(len(brickid))
    oldid = -1 
    image = None
    start = 0
    for i in range(len(brickid) + 1):
        if not (i == len(brickid) or brickid[i] != oldid): continue
        sl = slice(start, i)

        if i != 0:
            ind = (x[sl], y[sl])
            l = numpy.ravel_multi_index(ind, image.shape, mode='wrap')
            pixels[sl] = image.flat[l]
        
        # advance
        if i != len(brickid):
            try:
                image = fits.open(repo % dict(brickid=brickid[i]))[0].data
            except Exception as e:
                image = numpy.empty((1, 1))
                image[0, 0] = default
            oldid = brickid[i]
            start = i
    return pixels

def optimize(bid, ra, dec):
    arg = bid.argsort()
    invarg = numpy.empty_like(arg)
    invarg[arg] = numpy.arange(len(arg), dtype='i8')
    return ra[arg], dec[arg], invarg

def test398599():
    """ needs file coadd/image-398599-z.fits """
    bricks = fits.open('bricks.fits')
    bi = BrickIndex(bricks[1].data) 
    print 'testing on brick 398599'
    x, y = numpy.indices((3600, 3600))
    x = numpy.ravel(x) + 0.5
    y = numpy.ravel(y) + 0.5
    ra, dec = bi.revert([398599] * len(x), x, y)
    bid = bi.query_brick(ra, dec)
    ra, dec, invarg = optimize(bid, ra, dec)
    print 'unique bid', len(numpy.unique(bid))
    bxy = bi.query(ra, dec)
    img = load('coadd/image-%(brickid)d-z.fits', *bxy)
    print (~numpy.isnan(img)).sum()
    img2 = fits.open('coadd/image-398599-z.fits')[0].data[:]
    diff = img[..., invarg].reshape(3600, 3600) - img2
    assert (diff[400:-400, 400:-400] == 0).all()
    print 'passed'

if __name__ == '__main__':
    bricks = fits.open('bricks.fits')
    bi = BrickIndex(bricks[1].data) 
    print bricks[1].data[398599 - 1]
    test398599()
    #print bricks[1].data.dtype
    #print bricks[1].data[900]

#    print load('coadd/depth-%(brickid)d-z.fits.gz', *bxy)

    dec = (numpy.random.random(size=20000) - 0.5)* 10 + 10.
    ra = numpy.random.random(size=20000) * 360. 
    bid = bi.query_brick(ra, dec)
    print 'unique bid', len(numpy.unique(bid))
    ra, dec, invarg = optimize(bid, ra, dec)
    bxy = bi.query(ra, dec)
    print len(bxy.T)
    print numpy.isnan(load('coadd/depth-%(brickid)d-z.fits.gz', *bxy)).sum()
    #bi.test()

