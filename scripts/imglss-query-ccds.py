
from __future__ import print_function
import numpy as np
import h5py

from imaginglss             import DECALS
from imaginglss.model             import dataproduct
from imaginglss.analysis    import completeness
from imaginglss.utils       import output
from imaginglss.model.ccdtable import CCDTable

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

cli.add_argument("ccdattrs", type=lambda x: x.upper(), nargs="+",
        help="column name to query ")

ns = cli.parse_args()
decals = DECALS(ns.conf)

np.seterr(divide='ignore', invalid='ignore')

def main():
    ccdtable = CCDTable(ns.ccdfile)

    ccdtree = points(ccdtable.RA, ccdtable.DEC).tree
    with h5py.File(ns.query) as ff:
        RA = ff['RA'][:]
        DEC = ff['DEC'][:]

    print("size of query is %d" % len(RA))
    print("number of CCDS is %d" % len(ccdtable))
    attrnames = list(ccdtable.data.dtype.names) + ['NEXP','NEXP_G','NEXP_R','NEXP_Z']
    if len(ns.ccdattrs) == 0 or \
       any(attr not in attrnames for attr in ns.ccdattrs):
        raise RuntimeError("ccdattr not found, available ones are %s"
                    % str(attrnames))

    attrdtype = np.dtype([
        (attr, 'f8')
            for attr in ns.ccdattrs])

    r1 = np.zeros_like(RA, dtype=attrdtype)
    r2 = np.zeros_like(RA, dtype=attrdtype)
    N  = np.zeros_like(RA, dtype='f8')
    NG = np.zeros_like(RA, dtype='f8')
    NR = np.zeros_like(RA, dtype='f8')
    NZ = np.zeros_like(RA, dtype='f8')

    querytree = points(RA, DEC).tree
    
    def process(r, i, j):
        mask = ccdtable.query_inside(i, RA[j], DEC[j])
        i = i[mask]
        j = j[mask]

        v = ccdtable.data[i]
        for attr in ns.ccdattrs:
            if attr not in ['NEXP','NEXP_G','NEXP_R','NEXP_Z']:
                np.add.at(r1[attr], j, v[attr])
                np.add.at(r2[attr], j, v[attr] ** 2)
        np.add.at(N, j, 1)
        np.add.at(NG, j[v['FILTER'] == bytes('g','utf-8')], 1)
        np.add.at(NR, j[v['FILTER'] == bytes('r','utf-8')], 1)
        np.add.at(NZ, j[v['FILTER'] == bytes('z','utf-8')], 1)

    ccdtree.root.enum(querytree.root, np.radians(0.2), process)

    for attr in ns.ccdattrs:
        if attr not in ['NEXP','NEXP_G','NEXP_R','NEXP_Z']:
            r1[attr] /= N
        if attr == 'NEXP':
            r1[attr] = N
        if attr == 'NEXP_G':
            r1[attr] = NG
        if attr == 'NEXP_R':
            r1[attr] = NR
        if attr == 'NEXP_Z':
            r1[attr] = NZ
            
        if ((NG + NR + NZ) != N).sum() != 0:
            print('Exposures do not add up...')

        print('mean and std of %s from the query is %g %g' % (attr, np.nanmean(r1[attr]), np.nanstd(r1[attr])))

    with h5py.File(ns.query, 'r+') as ff:
        for attrname in ns.ccdattrs:
            COLUMNNAME = 'CCD-%s' % attrname
            if COLUMNNAME in ff:
                del ff[COLUMNNAME]
            ds = ff.create_dataset(COLUMNNAME, data=r1[attrname])
            ds.attrs.update(cli.prune_namespace(ns))

            print('written %s as %s in %s' % (attrname, COLUMNNAME, ns.query))

if __name__ == "__main__":
    main()
