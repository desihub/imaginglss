import numpy
from ..utils import fits

__author__ = "Yu Feng and Martin White"
__version__ = "1.0"
__email__  = "yfeng1@berkeley.edu or mjwhite@lbl.gov"


class Tycho(numpy.ndarray):
    """ Representing a Tycho catalogue that
        is used to veto objects near stars 
    """
    def __new__(kls, path):
        data = fits.read_table(path).view(type=kls)
        self = numpy.empty(len(data),
            dtype=[
                ('RA', 'f8'),
                ('DEC', 'f8'),
                ('VTMAG', 'f8'),
                ('VMAG', 'f8'),
                ('BMAG', 'f8'),
                ('BTMAG', 'f8'),
                ('VARFLAG', 'i8'),
                ])
        self['RA'] = data['RA']
        self['DEC'] = data['DEC']
        self['VARFLAG'] = data['VARFLAG']
        self['VMAG'] = data['VMAG']
        self['BMAG'] = data['BMAG']
        v = self['VMAG'] 
        b = self['BMAG'] 
        vt = v + 0.09 / 0.85 * (b - v)
        bt = 1.09 / 0.85 * (b - v) + v
        self['VTMAG'] = vt
        self['BTMAG'] = bt

        return self.view(type=kls)

