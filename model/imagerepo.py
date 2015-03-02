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

    def open(self, brick, **kwargs):
        kwargs['brickid'] = brick.id
        kwargs['brickname'] = brick.name
        return fits.open(self.get_filename(
            **kwargs
            ))[0].data[:]

    def get_filename(self, **kwargs):
        return os.path.join(self.root, 
            self.pattern) % kwargs

