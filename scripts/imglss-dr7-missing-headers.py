from fitsio import FITS
import os.path
import numpy
from argparse import ArgumentParser
from time import time

def get_meta(row, prefix):
    filename = row['IMAGE_FILENAME'].decode().strip()
    skycounts = row['CCDSKYCOUNTS']
    pixscale = row['PIXSCALE_MEAN']
    exptime = row['EXPTIME']
    ccdzpt = row['CCDZPT']
    fullname = os.path.join(prefix, filename)

    d = {}
    with FITS(fullname, upper=True) as f:
        h = f[0].read_header()
        d['SEEING'] = float(h['DIMMSEE'])
        d['AIRMASS'] = float(h['AIRMASS'])
        d['CCDSKYMAG'] = -2.5 * numpy.log10(skycounts / pixscale ** 2 / exptime) + ccdzpt

    return d

import sharedmem
import argparse
ap = argparse.ArgumentParser()

ap.add_argument('outputprefix', help='output file prefix. will write to ccds-annotated-dr7-missing.fits.gz')
ap.add_argument('--limit', help='limit to processing only N ccds', type=int, default=None)
ap.add_argument('--np', help='number of processes', type=int, default=8)
ap.add_argument('--prefix', help='prefix to dr7', default='/global/project/projectdirs/cosmo/data/legacysurvey/dr7')
ap.add_argument('--stagingprefix', help='prefix to staging files (ccds)',
            default='/global/project/projectdirs/cosmo/staging')

ns = ap.parse_args()

def main(ns):

    with FITS(os.path.join(ns.prefix, 'ccds-annotated-dr7.fits.gz'), upper=True) as f:

        print(f[1]['DEC'][:].ndim)
        print(f[1]['DEC'][3].ndim)

        if ns.limit is None:
            ann = f[1][:]
        else:
            ann = f[1][:ns.limit]

    result = sharedmem.empty(len(ann), dtype=
                [
                    ('SEEING', 'f8'),
                    ('AIRMASS', 'f8'),
                    ('CCDSKYMAG', 'f8'),
                ])
    Ntotal = [0]
    time_start = time()

    print('Number of ccd files', len(ann))
    with sharedmem.MapReduce(np=ns.np) as pool:
        chunksize = 128

        def work(i):
            for result_row, row in zip(result[i:i+chunksize], ann[i:i+chunksize]):
                d = get_meta(row, ns.stagingprefix)
                for key in d:
                    result_row[key] = d[key]
            return len(result[i:i+chunksize])

        def reduce(n):
            Ntotal[0] = Ntotal[0] + n
            usedtime = time() - time_start
            rate = Ntotal[0]  / usedtime
            eta = (len(ann) - Ntotal[0]) / rate
            print('%08d / %08d' % (Ntotal[0], len(ann)), 'ccds are scanned rate=%g eta=%g seconds' % (rate, eta))

        pool.map(work, range(0, len(ann), chunksize), reduce=reduce)

    with FITS(os.path.join(ns.outputprefix, 'ccds-annotated-missing-dr7.fits.gz'), upper=True, clobber=True, mode='rw') as f:
        f.write_table(result)


main(ns)

