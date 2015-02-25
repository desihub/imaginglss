import os
import os.path
import decamquery as dq
from astropy.io import fits
import numpy

class DECALS(object):
    def __init__(self, root=None):
        if root is None:
            root = os.env.get("DECALS_IMAGING", '.') 
        self.root = root
        bricks = fits.open(os.path.join(self.root, 'bricks.fits'))
        self.brickindex = dq.BrickIndex(bricks[1].data)

        self.bands = ['g','r','z']

        self.repos = dict(
            DEPTH='coadd/depth-%%(brickid)d-%(band)s.fits.gz',
            IMAGE='coadd/image-%%(brickid)d-%(band)s.fits',
            MODEL='coadd/model-%%(brickid)d-%(band)s.fits',
        )

        catalogue = numpy.array(fits.open(os.path.join(self.root, 'tractor-edr.fits'))[1].data, 
            copy=True)
        mask = catalogue['RA'] != 360
        self.catalogue = catalogue[mask]

    def get_pixel(self, RA, DEC, *args):
        RA, DEC, invarg = self.brickindex.optimize(RA, DEC)
        bxy = self.brickindex.query(RA, DEC)

        images = numpy.empty((len(RA), len(args)))
        for (i, (repo, band)) in enumerate(args):
            repo = self.repos[repo] % dict(band=band)
            repo = os.path.join(self.root, repo)
            img = dq.load(repo, *bxy)
            images[:,  i] = img[invarg]

        return images

if __name__ == '__main__':
    decals = DECALS('.')
    print len(decals.catalogue)
    dec = decals.catalogue['DEC'][:10000:10]
    ra = decals.catalogue['RA'][:10000:10]
    print decals.get_pixel(ra, dec, 
        ('DEPTH', 'g'),
        ('DEPTH', 'r'),
        ('DEPTH', 'z'),
        )
