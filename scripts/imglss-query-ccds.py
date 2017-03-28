
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
        help="survey ccd fits file")

cli.add_argument("ccdattr",
        help="column name to query ")

ns = cli.parse_args()
decals = DECALS(ns.conf)

np.seterr(divide='ignore', invalid='ignore')

import fitsio

class CCDTable(object):
    def __init__(self, filename):
        ccdhdu = fitsio.FITS(ns.ccdfile)[1]
        self.RA = ccdhdu['RA'][:]
        self.DEC = ccdhdu['DEC'][:]
        self.CD = np.array([
                ccdhdu['CD1_1'][:],
                ccdhdu['CD1_2'][:],
                ccdhdu['CD2_1'][:],
                ccdhdu['CD2_2'][:],]).copy()
        self.CRVAL = np.array([
                ccdhdu['CRVAL1'][:],
                ccdhdu['CRVAL2'][:]]).copy()
        # offset by one since Python starts from 0.
        self.CRPIX = np.array([
                ccdhdu['CRPIX1'][:],
                ccdhdu['CRPIX2'][:]]).copy() - 1
        self.SIZE = np.array([
                ccdhdu['WIDTH'][:],
                ccdhdu['HEIGHT'][:]]).copy()
        self.data = ccdhdu[:]

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

    if ns.ccdattr not in ccdtable.data.dtype.names:
        raise RuntimeError("ccdattr not found, available ones are %s"
                    % str(ccdtable.data.dtype.names))
        
    r1 = np.zeros_like(RA)
    r2 = np.zeros_like(RA)
    N = np.zeros_like(RA)

    querytree = points(RA, DEC).tree
    
    def process(r, i, j):
        mask = ccdtable.query_inside(i, RA[j], DEC[j])
        i = i[mask]
        j = j[mask]

        v = ccdtable.data[ns.ccdattr][i]
        np.add.at(r1, j, v)
        np.add.at(r2, j, v ** 2)
        np.add.at(N, j, 1)

    ccdtree.root.enum(querytree.root, np.radians(0.2), process)

    COLUMNNAME = 'CCD-%s' % ns.ccdattr.upper()

    r1 /= N

    print('mean and std of %s from the query is %g %g' % (ns.ccdattr, r1.mean(), r1.std()))

    with h5py.File(ns.query, 'r+') as ff:
        if COLUMNNAME in ff:
            del ff[COLUMNNAME]
        ds = ff.create_dataset(COLUMNNAME, data=r1)
        ds.attrs.update(cli.prune_namespace(ns))

    print('written as %s in %s' % ( COLUMNNAME, ns.query))

if __name__ == "__main__":
    main()
