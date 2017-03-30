
from __future__ import print_function
import numpy as np
import h5py

from imaginglss             import DECALS
from imaginglss.model             import dataproduct
from imaginglss.analysis    import completeness
from imaginglss.utils       import output
from imaginglss.utils.wcs_tangent import ang2pix

from imaginglss.cli import CLI
from kdcount.sphere import points

cli = CLI("Query completeness",
        enable_target_plugins=True,
        enable_confidence=True,
        enable_tycho_veto=True)

cli.add_argument("query",
        help="catalogue to query ccd systematic")

cli.add_argument("ccdfile",
        help="survey ccd fits file; this is the path to the -decals file. we also need the -nondecals and -extras file.")

cli.add_argument("ccdattr", type=lambda x: x.upper(),
        help="column name to query ")

ns = cli.parse_args()
decals = DECALS(ns.conf)

np.seterr(divide='ignore', invalid='ignore')

import fitsio

class CCDTable(object):
    def __init__(self, filename):
        data = []
        for fn in [filename,
                   filename.replace('ccds-decals', 'ccds-nondecals'),
                   filename.replace('ccds-decals', 'ccds-extra')]:
            ccdhdu = fitsio.FITS(fn, upper=True)[1]
            data.append(ccdhdu[:])

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

def main():
    ccdtable = CCDTable(ns.ccdfile)

    ccdtree = points(ccdtable.RA, ccdtable.DEC).tree
    with h5py.File(ns.query) as ff:
        RA = ff['RA'][:]
        DEC = ff['DEC'][:]

    print("size of query is %d" % len(RA))
    print("number of CCDS is %d" % len(ccdtable))

    if ns.ccdattr not in ccdtable.data.dtype.names and ns.ccdattr != 'NEXP':
        raise RuntimeError("ccdattr not found, available ones are %s"
                    % str(list(ccdtable.data.dtype.names) + ['NEXP']))
        
    r1 = np.zeros_like(RA)
    r2 = np.zeros_like(RA)
    N = np.zeros_like(RA)

    querytree = points(RA, DEC).tree
    
    def process(r, i, j):
        mask = ccdtable.query_inside(i, RA[j], DEC[j])
        i = i[mask]
        j = j[mask]

        if ns.ccdattr != 'NEXP':
            v = ccdtable.data[ns.ccdattr][i]
            np.add.at(r1, j, v)
            np.add.at(r2, j, v ** 2)
        np.add.at(N, j, 1)

    ccdtree.root.enum(querytree.root, np.radians(0.2), process)

    COLUMNNAME = 'CCD-%s' % ns.ccdattr

    if ns.ccdattr != 'NEXP':
        r1 /= N
    else:
        r1 = N

    print('mean and std of %s from the query is %g %g' % (ns.ccdattr, r1.mean(), r1.std()))

    with h5py.File(ns.query, 'r+') as ff:
        if COLUMNNAME in ff:
            del ff[COLUMNNAME]
        ds = ff.create_dataset(COLUMNNAME, data=r1)
        ds.attrs.update(cli.prune_namespace(ns))

    print('written as %s in %s' % ( COLUMNNAME, ns.query))

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

if __name__ == "__main__":
    main()
