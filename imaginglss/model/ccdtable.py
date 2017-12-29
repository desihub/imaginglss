from __future__ import print_function
import numpy as np

from ..utils.wcs_tangent import ang2pix
from ..utils import fits

class CCDTable(object):
    def __init__(self, filepath, filenames):
        data = []
        for fnrow in filenames:
            if not isinstance(fnrow, (tuple, list)):
                fnrow = (fnrow,)
            datarow = []
            for fn in fnrow:
               fn = filepath + fn
               try:
                   ccdhdu = fits.read_table(fn)
               except Exception as e:
                   print("CCDTable: file %s skipped due to %s" % (fn, e))
               datarow.append(ccdhud)

            data.append(datarow)

        data = concatenate_struct_arrays(data)

        mask = get_legacy_cataloged(data)

        self.data = data[mask]

        self.RA = data['RA'][mask]
        self.DEC = data['DEC'][mask]
        self.CD = np.array([
                data['CD1_1'][mask],
                data['CD1_2'][mask],
                data['CD2_1'][mask],
                data['CD2_2'][mask],]).copy()
        self.CRVAL = np.array([
                data['CRVAL1'][mask],
                data['CRVAL2'][mask]]).copy()
        # offset by one since Python starts from 0.
        self.CRPIX = np.array([
                data['CRPIX1'][mask],
                data['CRPIX2'][mask]]) - 1
        self.SIZE = np.array([
                data['WIDTH'][mask],
                data['HEIGHT'][mask]]).copy()

    def __len__(self):
        return len(self.data)

    def query_inside(self, ccdid, RA, DEC):
        CD = self.CD[..., ccdid]
        CRVAL = self.CRVAL[..., ccdid]
        CRPIX = self.CRPIX[..., ccdid]
        SIZE = self.SIZE[..., ccdid]
        xy = ang2pix((RA, DEC), CD, CRPIX, CRVAL)
        inside = ((xy >= 0) & (xy < SIZE)).all(axis=0)
        return inside

def concatenate_struct_arrays(arrays):
    """ sorry, but array shall be a list of list of arrays. (nested) """
    names = None
    for arow in arrays:
        for array in arow:
            if names is None:
                names = set(array.dtype.names)
            else:
                names = names.union(set(array.dtype.names))
    dtype = []
    for name in names:
        # only check the first row
        for array in arrays[0]:
            if name not in array.dtype[name]: continue
            dtype.append((name, array[0].dtype[name]))
            break
    dtype = np.dtype(dtype)
    result = np.empty(sum(len(arow[0]) for arow in arrays), dtype=dtype)
    offset = 0
    for arow in arrays:
        # use len of the first array in the nested row
        alen = len(arow[0])
        for name in names:
            for array in arow:
                if name not in array.dtype[name]: continue
                result[offset:offset+alen][name] = array[name]
                break
        offset = offset + alen

    return result

def get_legacy_cataloged(data):

    zp0 = dict(g=25.08, r=25.29, z=24.92)
    
    mask = data['EXPTIME'] >= 30
    mask &= data['CCDNMATCH'] >= 20
    mask &= np.abs(data['ZPT'] - data['CCDZPT']) <= 0.1

    for band in 'grz':
        hasband = data['FILTER'] == bytes(band,'utf-8')
        mask[hasband] &= data['ZPT'][hasband] >= (zp0[band] - 0.5) 
        mask[hasband] &= data['ZPT'][hasband] <= (zp0[band] + 0.25)

    kept = mask.sum() * 100.0 / len(mask)
    print('Percentage of CCDs propagated to Legacy catalogs:',np.round(kept,1))

    return mask
