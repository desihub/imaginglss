"""
This module provides an interface to query coadd images from 
DECAM imaging data release.

Will need bricks.fits and coadd/* (depending on which image to query)

    usage:

        bricks = fits.open('bricks.fits')
        bi = BrickIndex(bricks[1].data) 
        ra = ....
        dec = ....
        ra, dec, invarg = bi.optimize(ra, dec)
        brickid, x, y = bi.query(ra, dec)
        value = load('coadd/images-%(brickid)d-z.fits', brickid, x, y)

        # if we want value in original ra, dec order
        value = value[invarg]
 
"""
from astropy import wcs
from astropy.io import fits
import numpy

__all__ = [
    'BrickIndex',
    'load',
]

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

    def query_brick(self, RA, DEC):
        """ 
            finds BRICKID for given RA and DEC 

        """
        RA = numpy.asarray(RA)
        DEC = numpy.asarray(DEC)
        row = numpy.int32(numpy.floor(DEC * self.ROWMAX / 180 + 360. + 0.5))
        row = numpy.clip(row, 0, self.ROWMAX)
        ncols = self.ncols[row]
        col = numpy.int32(numpy.floor(RA * ncols / 360. ))
        hash = row * (self.COLMAX + 1) + col
        ind = self.hash.searchsorted(hash)
        bid = self.hdudata['BRICKID'][ind]
        return bid

    def optimize(self, ra, dec, active_bid=None):
        """ optimize the ordering of ra, dec,

            This will sort ra, dec by their brickid.

            return sorted_ra, sorted_dec, invert_arg

            invariance:
                sorted_ra[inverg_arg] == ra
                sorted_dec[inverg_arg] == dec
        """
        bid = self.query_brick(ra, dec)
        arg = bid.argsort()

        invarg = numpy.empty_like(arg)
        invarg[arg] = numpy.arange(len(arg), dtype='i8')
        return ra[arg], dec[arg], invarg


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
        pix[:, 0] = brk

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
                ra = self.hdudata['RA'][brk[i] - 1]  # watch out
                dec = self.hdudata['DEC'][brk[i] - 1] # watch out
                header['CRVAL1']  =     ra, # / Reference RA                                   
                header['CRVAL2']  =     dec, # / Reference Dec                                  

                q = wcs.WCS(header)

                oldbrk = brk[i]
                start = i

        return pix.T

def load(repo, brickid, x, y, default=numpy.nan):
    """ 
        Load all pixels indexed by brickid, x, y from the primary HDU of fits files in repo,

        repo is a format string with %(brickid).

        eg:
            bxy = bi.query([243.6] * 6, numpy.arange(11.75 - 0.12, 11.75 + 0.12, 0.04))
            print load('coadd/depth-%(brickid)d-z.fits.gz', *bxy.T)

        This is faster if brickid, x, y is sorted by brickid.

        Note that brickid starts from 1.

        Of course the files has to be there already!
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
                f = fits.open(repo % dict(brickid=brickid[i]))
                print 'opening file', brickid[i]
                image = numpy.array(f[0].data, copy=True)
                f.close()
            except Exception as e:
                image = numpy.empty((1, 1))
                image[0, 0] = default
            oldid = brickid[i]
            start = i
    return pixels

def test398599():
    """ test image readout on brick-398599. 
        needs file coadd/image-398599-z.fits
        this test took ~21 seconds
    """
    bricks = fits.open('bricks.fits')
    bi = BrickIndex(bricks[1].data) 
    print 'testing on brick 398599'
    x, y = numpy.indices((3600, 3600))
    x = numpy.ravel(x) + 0.5
    y = numpy.ravel(y) + 0.5
    ra, dec = bi.revert([398599] * len(x), x, y)
    ra, dec, invarg = bi.optimize(ra, dec)
    bxy = bi.query(ra, dec)
    img = load('coadd/image-%(brickid)d-z.fits', *bxy)
    print (~numpy.isnan(img)).sum()
    img2 = fits.open('coadd/image-398599-z.fits')[0].data[:]
    img = img[..., invarg]
    diff = img.reshape(3600, 3600) - img2

    # FIXME: tighten this up
    assert (diff[300:-300, 300:-300] == 0).all()
    print 'passed'
def testquery_brick():
    print 'testing query_brick'
    bricks = fits.open('bricks.fits')
    bi = BrickIndex(bricks[1].data) 
    
    rows = bi.hdudata[bi.query_brick(bi.hdudata['RA'], bi.hdudata['DEC']) - 1]
    print 'id of failed:', (rows['BRICKID'] != bi.hdudata['BRICKID']).nonzero()
    rows = bi.hdudata[bi.query_brick(bi.hdudata['RA'], bi.hdudata['DEC'] - 0.125) - 1]
    print 'id of failed:', (rows['BRICKID'] != bi.hdudata['BRICKID']).nonzero()
    rows = bi.hdudata[bi.query_brick(bi.hdudata['RA'], bi.hdudata['DEC'] + 0.124) - 1]
    print 'id of failed:', (rows['BRICKID'] != bi.hdudata['BRICKID']).nonzero()
    print 'assert seeing [] [] [] above'

if __name__ == '__main__':
    testquery_brick()
    test398599()
    if False:
        bricks = fits.open('bricks.fits')
        bi = BrickIndex(bricks[1].data) 
        print bricks[1].data[398599 - 1]
        dec = (numpy.random.random(size=20000) - 0.5)* 10 + 10.
        ra = numpy.random.random(size=20000) * 360. 
        ra, dec, invarg = bi.optimize(ra, dec)
        bxy = bi.query(ra, dec)
        print len(bxy.T)
        print numpy.isnan(load('coadd/depth-%(brickid)d-z.fits.gz', *bxy)).sum()
    #bi.test()

