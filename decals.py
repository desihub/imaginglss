import os
import os.path
import decamquery as dq
from astropy.io import fits
import numpy
import glob
import re
def contains(haystack, needle):
    """ returns true if needle is in haystack """

    ind = haystack.searchsorted(needle)
    ind.clip(0, len(haystack) - 1, ind)
    return haystack[ind] == needle

class DECALS(object):
    def __init__(self, root=None):
        if root is None:
            root = os.env.get("DECALS_IMAGING", '.') 
        self.root = root

        bricks = fits.open(os.path.join(self.root, 'bricks.fits'))[1].data
        self.brickindex = dq.BrickIndex(bricks)

        self.bands = {'u':0, 'g':1,'r':2,'i':3, 'z':4, 'Y':5}

        self.repos = dict(
            DEPTH='coadd/depth-%%(brickid)d-%(band)s.fits.gz',
            IMAGE='coadd/image-%%(brickid)d-%(band)s.fits',
            MODEL='coadd/model-%%(brickid)d-%(band)s.fits',
        )

        observed_bricks = numpy.unique([
            int(re.search('-([0123456789]+)\.', fn).group(1))
                for fn in glob.glob(os.path.join(self.root, 'tractor/tractor-[0-9]*.fits'))
        ])

        catalogue = []
        for brick in observed_bricks:
            catalogue.append(numpy.array(
                fits.open(os.path.join(self.root, 'tractor/tractor-%d.fits' % brick))[1].data, 
                    copy=True))
        catalogue = numpy.concatenate(catalogue)
        #mask = catalogue['RA'] != 360
        catalogue['RA'] %= 360.

        self.catalogue = catalogue #[mask]
        self.galaxies = catalogue[catalogue['TYPE'] != 'S']


        self.observed_bricks = observed_bricks

        # approximate area in degrees. Currently a brick is 0.25 * 0.25 deg**2
        self.observed_area = 41253. * len(self.observed_bricks) / len(bricks)

    def get_pixel(self, RA, DEC, keys, default=numpy.nan):
        images = numpy.empty((len(RA), len(keys)))

        images[...] = default

        bid = self.brickindex.query_brick(RA, DEC)
        mask = contains(self.observed_bricks, bid)
        ra = RA[mask]
        dec = DEC[mask]
        ra, dec, invarg = self.brickindex.optimize(ra, dec)
        bxy = self.brickindex.query(ra, dec)

        print 'querying', len(ra), 'out of', len(RA)
        for (i, (repo, band)) in enumerate(keys):
            repo = self.repos[repo] % dict(band=band)
            repo = os.path.join(self.root, repo)
            img = dq.load(repo, *bxy, default=default)
            images[:,  i][mask] = img[invarg]

        return images

if __name__ == '__main__':
    decals = DECALS('.')
    print len(decals.catalogue)
    
    dec = decals.catalogue['DEC']#[:1000000:]
    ra = decals.catalogue['RA']#[:100000:]
    image = decals.get_pixel(ra, dec, 
        [('DEPTH', 'g'),
        ('DEPTH', 'r'),
        ('DEPTH', 'z'),]
        )
    numpy.savetxt('depth.txt', image) 
