import os
import os.path
from astropy.io import fits
import numpy
import glob
import re

from . import brickindex
from . import imagerepo
class Lazy(object):
    def __init__(self, calculate_function):
        self._calculate = calculate_function

    def __get__(self, obj, _=None):
        if obj is None:
            return self
        value = self._calculate(obj)
        setattr(obj, self._calculate.func_name, value)
        return value
    
def contains(haystack, needle):
    """ returns true if needle is in haystack """

    ind = haystack.searchsorted(needle)
    ind.clip(0, len(haystack) - 1, ind)
    return haystack[ind] == needle

class DataRelease(object):
    def __init__(self, root=None, version=None):
        if root is None:
            root = os.env.get("DECALS_IMAGING", '.') 
        self.root = root

        bricks = fits.open(os.path.join(self.root, 'bricks.fits'))[1].data
        self.brickindex = brickindex.BrickIndex(bricks)

        self.bands = {'u':0, 'g':1, 'r':2, 'i':3, 'z':4, 'Y':5}

        self.images = dict(
            DEPTH=imagerepo.ImageRepo(self.root, 'coadd/depth-%(brickid)d-%(band)s.fits.gz'),
            IMAGE=imagerepo.ImageRepo(self.root, 'coadd/image-%(brickid)d-%(band)s.fits'),
            MODEL=imagerepo.ImageRepo(self.root, 'coadd/model-%(brickid)d-%(band)s.fits'),
        )

        observed_bricks = numpy.unique([
            int(re.search('-([0123456789]+)\.', fn).group(1))
                for fn in glob.glob(os.path.join(self.root, 'tractor/tractor-[0-9]*.fits'))
        ])
        self.observed_bricks = observed_bricks

        # approximate area in degrees. Currently a brick is 0.25 * 0.25 deg**2
        self.observed_area = 41253. * len(self.observed_bricks) / len(bricks)

    @Lazy
    def catalogue(self):
        catalogue = []
        for brick in observed_bricks:
            catalogue.append(numpy.array(
                fits.open(os.path.join(self.root, 'tractor/tractor-%d.fits' % brick))[1].data, 
                    copy=True))
        catalogue = numpy.concatenate(catalogue)
        catalogue['RA'] %= 360.

        return catalogue

    def readout(self, coord, keys, default=numpy.nan):
        """ readout pixels at coord.
            querying from several images. 
            
            example:
                corrd = (RA, DEC)
                keys = [('DEPTH', 'z'), ....]
        """
        RA, DEC = coord
        images = numpy.empty((len(RA), len(keys)))
        images[...] = default

        bid = self.brickindex.query((RA, DEC))
        mask = contains(self.observed_bricks, bid)
        ra = RA[mask]
        dec = DEC[mask]
        coord, invarg = self.brickindex.optimize((ra, dec))
        bid = self.brickindex.query(coord)

        pixels = numpy.empty(len(bid), 'f8')
        pixels[:] = default

        
        ubid = numpy.unique(bid)
        print ubid
        for (i, (repo, band)) in enumerate(keys):
            for b in ubid:
                brick = self.brickindex[b]
                first = bid.searchsorted(b, side='left')
                last = bid.searchsorted(b, side='right')
                sl = slice(first, last)

                x, y = numpy.int32(brick.query(coord[:, sl]))

                img = self.images[repo].open(brick, band=band)
                l = numpy.ravel_multi_index((x, y), img.shape, mode='raise')
                pixels[sl] = img.flat[l] 
            #
            images[:,  i][mask] = pixels[invarg]
            
        return images
