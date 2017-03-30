from __future__ import print_function
import numpy as np

from ..utils.wcs_tangent import ang2pix
from ..utils import fits

class CCDTable(object):
    def __init__(self, filename):
        data = []
        for fn in [filename,
                   filename.replace('ccds-decals', 'ccds-nondecals'),
                   filename.replace('ccds-decals', 'ccds-extra')]:
            try:
                ccdhdu = fits.read_table(fn)
            except Exception as e:
                print("CCDTable: file %s skipped due to %s" % (fn, e))

            data.append(ccdhdu)

        data = concatenate_struct_arrays(data)

        self.data = data

        self.RA = data['RA'][:]
        self.DEC = data['DEC'][:]
        self.CD = np.array([
                data['CD1_1'][:],
                data['CD1_2'][:],
                data['CD2_1'][:],
                data['CD2_2'][:],]).copy()
        self.CRVAL = np.array([
                data['CRVAL1'][:],
                data['CRVAL2'][:]]).copy()
        # offset by one since Python starts from 0.
        self.CRPIX = np.array([
                data['CRPIX1'][:],
                data['CRPIX2'][:]]) - 1
        self.SIZE = np.array([
                data['WIDTH'][:],
                data['HEIGHT'][:]]).copy()

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

