from __future__ import print_function
import numpy as np

from ..utils.wcs_tangent import ang2pix
from ..utils import fits

class CCDTable(object):
    def __init__(self, filepath, filenames):
        data = []
        for fn in filenames:
            fn = filepath + fn
            try:
                ccdhdu = fits.read_table(fn)
            except Exception as e:
                print("CCDTable: file %s skipped due to %s" % (fn, e))

            data.append(ccdhdu)

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
    names = None
    for array in arrays:
        if names is None:
            names = set(array.dtype.names)
        else:
            names = names.union(set(array.dtype.names))
    dtype = []
    for name in names:
        dtype.append((name, array[0].dtype[name]))
    dtype = np.dtype(dtype)
    result = np.empty(sum(len(a) for a in arrays), dtype=dtype)
    offset = 0
    for array in arrays:
        for name in names:
            result[offset:offset+len(array)][name] = array[name]
        offset = offset + len(array)
        
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
