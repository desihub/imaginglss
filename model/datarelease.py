import os
import os.path
from astropy.io import fits
import numpy
import glob
import re

from . import brickindex
from . import imagerepo
from . import catalogue

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
    def __init__(self, root=None, cacheroot=None, version=None):
        if root is None:
            root = os.environ.get("DECALS_IMAGING", '.') 
        self.root = root

        if cacheroot is None:
            cacheroot = os.environ.get("DECALS_CACHE", '.') 

        self.cacheroot = cacheroot

        bricks = fits.open(os.path.join(self.root, 'bricks.fits'))[1].data
        self.brickindex = brickindex.BrickIndex(bricks)

        self.bands = {'u':0, 'g':1, 'r':2, 'i':3, 'z':4, 'Y':5}

        self.images = {}
        for band in self.bands:
            self.images[band] = dict(
                DEPTH=imagerepo.ImageRepo(self.root, 
                    'coadd/depth-%%(brickid)d-%(band)s.fits.gz' % dict(band=band)),
                IMAGE=imagerepo.ImageRepo(self.root, 
                    'coadd/depth-%%(brickid)d-%(band)s.fits.gz' % dict(band=band)),
                MODEL=imagerepo.ImageRepo(self.root, 
                    'coadd/depth-%%(brickid)d-%(band)s.fits.gz' % dict(band=band)),
            )

        observed_brickids = numpy.unique([
            int(re.search('-([0123456789]+)\.', fn).group(1))
                for fn in glob.glob(os.path.join(self.root, 'tractor/tractor-[0-9]*.fits'))
        ])
        self.observed_brickids = observed_brickids

        # approximate area in degrees. Currently a brick is 0.25 * 0.25 deg**2
        self.observed_area = 41253. * len(self.observed_brickids) / len(bricks)

        self.catalogue = catalogue.Catalogue(
            os.path.join(self.cacheroot, 'catalogue'),
            [
            os.path.join(self.root, 'tractor/tractor-%d.fits' % brick)
            for brick in self.observed_brickids])
        # fix RA
        #self.catalogue['RA'][:] %= 360.
            
    def readout(self, coord, repos, default=numpy.nan):
        """ readout pixels at coord.
            querying from several images. 
            
            example:
                corrd = (RA, DEC)
                repos is fetched from self.images

            This is here, because we want to query multiple images
            at the same time. 
            It is also convenient to have it here to make use of
            brickindex. (ImageRepo is then just a stub with no business logic)

            Otherwise it makes more sense to
            have readout in ImageRepo.
        """
        RA, DEC = coord
        images = numpy.empty((len(RA), len(repos)))
        images[...] = default

        bid = self.brickindex.query((RA, DEC))
        # watch out bid + 1
        mask = contains(self.observed_brickids, bid + 1)
        ra = RA[mask]
        dec = DEC[mask]
        coord, invarg = self.brickindex.optimize((ra, dec))
        bid = self.brickindex.query(coord)

        pixels = numpy.empty(len(bid), 'f8')
        pixels[:] = default

        
        ubid = numpy.unique(bid)

        for (i, repo) in enumerate(repos):
            for b in ubid:
                brick = self.brickindex[b]
                first = bid.searchsorted(b, side='left')
                last = bid.searchsorted(b, side='right')
                sl = slice(first, last)

                x, y = numpy.int32(brick.query(coord[:, sl]))

                img = repo.open(brick)
                l = numpy.ravel_multi_index((x, y), img.shape, mode='raise')
                pixels[sl] = img.flat[l] 
            #
            images[:,  i][mask] = pixels[invarg]
            
        return images
