from fitsio import FITS
import os.path
import numpy
from argparse import ArgumentParser
from time import time

def get_meta(row, prefix, cache=None):
    filename = row['IMAGE_FILENAME'].decode().strip()
    if cache is not None and filename in cache:
        return cache[filename]

    skycounts = row['CCDSKYCOUNTS']
    hdu = row['IMAGE_HDU']
    pixscale = row['PIXSCALE_MEAN']
    exptime = row['EXPTIME']
    ccdzpt = row['CCDZPT']
    filter = row['FILTER'].decode().upper().strip()
    fullname = os.path.join(prefix, filename)

    d = {}
    with FITS(fullname, upper=True) as f:
        h = f[0].read_header()

        seeing = numpy.nan
        if 'DIMMSEE' in h:
            seeing = float(h['DIMMSEE'])
        if numpy.isnan(seeing) and 'DIMM2SEE' in h:
            seeing = float(h['DIMM2SEE'])

        d['SEEING'] = seeing

        if numpy.isnan(d['SEEING']):
            print(fullname, 'DIMMSEE is NaN')

        airmass = numpy.nan
        if 'AIRMASS' in h:
            airmass = float(h['AIRMASS'])

        d['AIRMASS'] = airmass

        if numpy.isnan(d['AIRMASS']):
            print(fullname, 'AIRMASS is nan')
        
        d['CCDSKYMAG'] = -2.5 * numpy.log10(skycounts / pixscale ** 2 / exptime) + ccdzpt

#        print(fullname, hdu, list(h.keys()))
    if cache is not None:
        cache[filename] = d
    return d

import sharedmem
import argparse
ap = argparse.ArgumentParser()

ap.add_argument('outputprefix', help='output file prefix. will write to ccds-annotated-dr7-missing.fits.gz')
ap.add_argument('--limit', help='limit to processing only N ccds', type=int, default=None)
ap.add_argument('--np', help='number of processes', type=int, default=8)
ap.add_argument('--prefix', help='prefix to dr7', default='/global/project/projectdirs/cosmo/data/legacysurvey/dr7')
ap.add_argument('--dr5-prefix', help='prefix to dr5', default='/global/project/projectdirs/cosmo/data/legacysurvey/dr5')
ap.add_argument('--stagingprefix', help='prefix to staging files (ccds)',
            default='/global/project/projectdirs/cosmo/staging')

ns = ap.parse_args()

def main(ns):

    with FITS(os.path.join(ns.prefix, 'ccds-annotated-dr7.fits.gz'), upper=True) as f:
        if ns.limit is None:
            ann = f[1][:]
        else:
            ann = f[1][:ns.limit]

        print('unique dr7 filenames', len(numpy.unique(ann['IMAGE_FILENAME'])))
        CCDNAMES = numpy.unique(ann['CCDNAME'])
        CCDNAMES.sort()
        FILTERS = numpy.unique(ann['FILTER'])
        FILTERS.sort()
        ann_keys = (ann['EXPNUM'] * 64 + CCDNAMES.searchsorted(ann['CCDNAME'])) * 8 + FILTERS.searchsorted(ann['FILTER'])

    with FITS(os.path.join(ns.dr5_prefix, 'ccds-annotated-dr5.fits.gz'), upper=True) as f:
        ann5 = f[1][:]

        print('unique dr5 filenames', len(numpy.unique(ann5['IMAGE_FILENAME'])))
        ann5_keys = (ann5['EXPNUM'].astype(numpy.int64) * 64 + CCDNAMES.searchsorted(ann5['CCDNAME'])) * 8 + FILTERS.searchsorted(ann5['FILTER'])
        arg = ann5_keys.argsort()
        ann5_keys = ann5_keys[arg]
        ann5 =  ann5[arg]

    print("dr7 has %d ccds" % len(ann))
    print("dr5 has %d ccds" % len(ann5))
    print(len(numpy.unique(ann_keys)))
    ind = ann5_keys.searchsorted(ann_keys)
    ind[ind < 0] = 0 
    ind[ind >= len(ann5_keys)] = 0 
    mask = ann5_keys[ind] == ann_keys
    print('found %d ccds' % mask.sum())
    assert len(numpy.unique(ann5['EXPNUM'] * 64 + ann5['CCDNUM'])) == len(ann5)
    #ann5.searchsorted(ann['EXPNUM'])
    result = sharedmem.empty(len(ann), dtype=
                [
                    ('SEEING', 'f8'),
                    ('AIRMASS', 'f8'),
                    ('CCDSKYMAG', 'f8'),
                ])
    Ntotal = [0]
    Nphys = [0]
    time_start = time()

    print('Number of ccd files', len(ann))
    cache = {}

    with sharedmem.MapReduce(np=ns.np) as pool:
        chunksize = 128

        def work(i):
            N = 0
            for result_row, row, key in zip(result[i:i+chunksize], ann[i:i+chunksize], ann_keys[i:i+chunksize]):
                ind = ann5_keys.searchsorted(key)
                if ind < 0 or ind >= len(ann5_keys): ind = 0

                if ann5_keys[ind] == key:
                    if (os.path.basename(ann5[ind]['IMAGE_FILENAME'].strip().decode()) 
                     == os.path.basename(row['IMAGE_FILENAME'].strip().decode())):
                        # found
                        for key in result.dtype.names:
                            result_row[key] = ann5[ind][key]
                        continue
    #                else:
    #                    print(ann5[ind]['IMAGE_FILENAME'].strip(), row['IMAGE_FILENAME'].strip())

                # get from file
                d = get_meta(row, ns.stagingprefix, cache)
                for key in d:
                    result_row[key] = d[key]
                N = N + 1
            return len(result[i:i+chunksize]), N

        def reduce(N, Nfile):
            Ntotal[0] = Ntotal[0] + N
            Nphys[0] = Nphys[0] + Nfile
            usedtime = time() - time_start
            rate = Ntotal[0]  / usedtime
            eta = (len(ann) - Ntotal[0]) / rate
            print('%08d / %08d' % (Ntotal[0], len(ann)), 'ccds are scanned rate=%g eta=%g seconds reused %d entries' % (rate, eta, Ntotal[0] - Nphys[0]))

        pool.map(work, range(0, len(ann), chunksize), reduce=reduce)

    with FITS(os.path.join(ns.outputprefix, 'ccds-annotated-missing-dr7.fits.gz'), upper=True, clobber=True, mode='rw') as f:
        f.write_table(result)


main(ns)

