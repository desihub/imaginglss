class Brick(object):
    def __init__(self, id, name, wcs):
        """ initialize a brick object 
            Brick is immutable .

            wcs describes the coordinate tranformation from 
            RA DEC to xy. 

            def wcs(ra, dec, xout, yout):
                ...
                xout[:], yout[:] = .....

        """
        self.id = id
        self.name = name
        self.wcs = wcs 

    def query(self, coord):
        """ returns the xy index of pixels for coord
            coord can be:
                (RA, DEC) tuple of arrays
                array of shape (2xN) RA DEC
            returns (2xN)
        """
        #FIXME: other types of input
        ra, dec = coord
        out = self.wcs.all_world2pix((ra, dec), 0).T
        return out
