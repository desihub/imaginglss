"""
This file generates the SFD map of E(B-V) from the standard
dust map DECALS uses.

It depends on astropy and fitsio explicitly.
model.datarelease dependency is used only to look up filenames
 
It depends on sfdmap.py, extracted from tractor source code.


usage:

make sure DECALS_IMAGING and DECALS_CACHING are initialized (see README.md)

python build_sfdmap.py
 Yu Feng <yfeng1@berkeley.edu>
 Martin White <mwhite@berkeley.edu>
"""

# TODO: skip already genearted files.
# TODO: use smaller pixels?
# FIXME: hard coded dustdir, only works on edison 

import fitsio
from model.sfdmap import SFDMap
import os.path
import os
import numpy
from model.utils import wcs_tangent
sfdmap = SFDMap(dustdir='/project/projectdirs/desi/software/edison/dust/v0_0/')

def process_file(header, newfilename=None):
    fac = 8
    shape = int(header['NAXIS1']) / fac, int(header['NAXIS2']) / fac
    header['NAXIS1'] = shape[0]
    header['NAXIS2'] = shape[1]
    header['CRPIX1'] = shape[0] / 2 + 0.5
    header['CRPIX2'] = shape[1] / 2 + 0.5
    header['CD1_1'] = float(header['CD1_1']) * fac
    header['CD2_2'] = float(header['CD2_2']) * fac
    header['IMTYPE'] = 'ebv'
    header.pop('FILTER', None)
    y, x = numpy.array(numpy.indices(shape), dtype='f8').reshape(2, -1)
    x += 0.5
    y += 0.5
    world = wcs_tangent.pix2ang_hdr(numpy.array((x,y)), header, zero_offset=True)
    ebv, junk = sfdmap.extinction([], world[0], world[1], get_ebv=True)
    ebv = ebv.reshape(shape)
    print newfilename, world.max(axis=1), world.min(axis=1)
    try:
        os.unlink(newfilename)
    except OSError:
        pass
    try:
        os.makedirs(os.path.dirname(newfilename))
    except OSError:
        pass
    fitsio.write(newfilename, ebv, header=header)

def main():
    from model.datarelease import DataRelease
    dr = DataRelease()
    for brick in dr.observed_bricks:
        header = dr.images['depth']['r'].metadata(brick)
        newfilename = dr.images['ebv'].get_filename(brick)
        if os.path.exists(newfilename):
            return
        process_file(header, newfilename)

def nomodel():
    root = os.environ['DECALS_IMAGING']
    for root, dirnames, filenames in \
        os.walk(os.path.join(root, 'coadd'), followlinks=True):
        for filename in filenames:
            if 'depth' not in filename: continue
            if '-r.fits' not in filename: continue
            filename = os.path.join(root, filename)
            header = dict(fitsio.read_header(filename))
            brickname = header['BRICKNAM']
            newfilename = 'decals-%s-ebv.fits' % brickname
            newfilename = os.path.join('dust', newfilename)
            process_file(header, newfilename)

from sys import argv
if len(argv) > 1 and argv[1] == '--no-model':
    nomodel()
else:
    from model.datarelease import DataRelease
    main()
