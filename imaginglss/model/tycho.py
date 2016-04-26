import numpy
from ..utils import fits

__author__ = "Yu Feng and Martin White"
__version__ = "1.0"
__email__  = "yfeng1@berkeley.edu or mjwhite@lbl.gov"


class Tycho1(numpy.ndarray):
    """ Representing a Tycho catalogue that
        is used to veto objects near stars 
    """
    def __new__(kls, path):
        data = fits.read_table(path).view(type=kls)
        if 'VMAG' not in data.dtype.names:
            raise ValueError("The file at `%s` is not a legacysurvey Tycho catalogue. Retrive from https://s3-us-west-1.amazonaws.com/imaginglss/tycho.tar.gz. (And unzip)" % path)
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

class Tycho2(numpy.ndarray):
    """ Representing a Tycho catalogue that
        is used to veto objects near stars 
    """
    def __new__(kls, path):
        data = fits.read_table(path).view(type=kls)
        if 'MAG_VT' not in data.dtype.names:
            raise ValueError("The file at `%s` is not a legacysurvey Tycho2 catalogue. Retrive from https://s3-us-west-1.amazonaws.com/imaginglss/tycho2.tar.gz. (And unzip)" % path)
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
        self['VTMAG'] = data['MAG_VT']
        self['BTMAG'] = data['MAG_BT']
        vt = self['VTMAG']
        bt = self['BTMAG']
        b = vt - 0.09 * (bt - vt)
        v = b - 0.85 * (bt - vt)
        self['VMAG'] = v
        self['BMAG'] = b

        return self.view(type=kls)

# Use tycho2 catalogue.
Tycho = Tycho2
