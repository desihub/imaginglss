from astropy import wcs
import numpy
class Brick(object):
    def __init__(self, id, name, ra, dec, ra1, ra2, dec1, dec2):
        """ initialize a brick object 
            Brick is immutable .

            wcs describes the coordinate tranformation from 
            RA DEC to xy. 
            def wcs(ra, dec, xout, yout):
                ...
                xout[:], yout[:] = .....  """
        self.id = id
        self.name = name
        self.ra = ra
        self.dec = dec
        self.ra1 = ra1
        self.ra2 = ra2
        self.dec1 = dec1
        self.dec2 = dec2

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
        header['CRVAL1']  =     ra, # / Reference RA                                   
        header['CRVAL2']  =     dec, # / Reference Dec                                  
        q = wcs.WCS(header)
        self.wcs = q

    def __repr__(self):
        return ("Brick(id=%d, name=%s, ra=%g, dec=%g, ...)"
            % (self.id, self.name, self.ra, self.dec ))

    def query(self, coord):
        """ returns the xy index of pixels for coord
            coord can be:
                (RA, DEC) tuple of arrays
                array of shape (2xN) RA DEC
            returns (2xN)
        """
        #FIXME: other types of input
        coord = numpy.array(coord).T
        out = self.wcs.all_world2pix(coord, 0).T
        return out
    def revert(self, xy):
        """ returns the RA, DEC index of pixels for coord
            coord can be:
                (x, y) tuple of arrays
                array of shape (2xN) x, y
            returns (2xN)
        """
        #FIXME: other types of input
        xy = numpy.array(xy).T
        out = self.wcs.all_pix2world(xy, 0).T
        return out
