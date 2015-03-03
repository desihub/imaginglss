from astropy.io import fits
import os.path

class ImageRepo(object):
    """ Image repository.
        Fetching image data from a brick

        pattern is "image-%(brick)s-r.fits.gz"
    """
    def __init__(self, root, pattern):
        self.root = root
        self.pattern = pattern
        self.cache = {}
        self.meta_cache = {} 

    def preload(self, bricks, **kwargs):
        for b in bricks:
            self.cache[b] = self.open(b, **kwargs)

    def open(self, brick, **kwargs):
        if brick in self.cache:
            return self.cache[brick]
        kwargs['brickid'] = brick.id
        kwargs['brickname'] = brick.name
        _ = self.metadata(brick, **kwargs)
        return fits.open(self.get_filename(
            **kwargs
            ))[0].data[:]

    def metadata(self, brick, **kwargs):
        if brick not in self.meta_cache:
            kwargs['brickid'] = brick.id
            kwargs['brickname'] = brick.name
            hdu = fits.open(self.get_filename(**kwargs))[0]
            meta = dict(hdu.header)
            self.meta_cache[brick] = meta
        return self.meta_cache[brick]

    def get_filename(self, **kwargs):
        return os.path.join(self.root, 
            self.pattern) % kwargs

