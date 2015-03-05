from astropy.io import fits
import os.path

class ImageRepo(object):
    """ Image repository.

        Image repository concerns about file system layout and
        IO of bricks. 

        Given a Brick object, ImageRepo can return the image data
        or meta data of the image.

    """
    def __init__(self, root, pattern):
        """
            pattern is a python formatting string.
            supported keywords are

            pattern can also be a callable with signature:

            def pattern(brick):
                return filename_without_root

            %(brickid)d
            %(brickname)d

            We need to extend this for DR1.
        """
        self.root = root
        self.pattern = pattern
        self.cache = {}
        self.meta_cache = {} 

    def preload(self, bricks, **kwargs):
        """ preload the images for 
            Brick objects listed in bricks list.

            **kwargs is unused 
        """
        for b in bricks:
            self.cache[b] = self.open(b, **kwargs) 
            
    def open(self, brick, **kwargs):
        """ open and read the image for 
            Brick object brick. The image content is returned
            as an array.

            **kwargs is unused 
        """
        if brick in self.cache:
            return self.cache[brick]
        _ = self.metadata(brick, **kwargs)
        return fits.open(self.get_filename(brick,
            **kwargs
            ))[0].data[:]

    def metadata(self, brick, **kwargs):
        """ Fetch the meta data about a brick.
        """
        if brick not in self.meta_cache:
            hdu = fits.open(self.get_filename(brick, **kwargs))[0]
            meta = dict(hdu.header)
            self.meta_cache[brick] = meta
        return self.meta_cache[brick]

    def get_filename(self, brick, **kwargs):
        if hasattr(self.pattern, '__call__'):
            return os.path.join(self.root, 
                self.pattern(brick, **kwargs))
        else:
            kwargs['brickid'] = brick.id
            kwargs['brickname'] = brick.name
            return os.path.join(self.root, 
                self.pattern) % kwargs

