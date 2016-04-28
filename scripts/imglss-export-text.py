import h5py
import numpy as np

from   imaginglss             import DECALS
from imaginglss.analysis import tycho_veto
from imaginglss.model import product

from argparse import ArgumentParser

ap = ArgumentParser()
ap.add_argument("catalogue", help="internal catalogue of HDF5 type.")
ap.add_argument("output", help="text file to store the catalogue." )
ap.add_argument("--conf", default=None,
        help="Path to the imaginglss config file, default is from DECALS_PY_CONFIG")

allvetos = [i for i in dir(tycho_veto) if not str(i).startswith( '_' )]
ap.add_argument("--use-tycho-veto", type=str, choices=allvetos, default=None)
ap.add_argument("--bands", nargs='+', type=str, choices=product.bands.keys(), default=['r', 'g', 'z'])
ap.add_argument("--sigma-z", type=float, default=3.0)
ap.add_argument("--sigma-g", type=float, default=5.0)
ap.add_argument("--sigma-r", type=float, default=5.0)

ns = ap.parse_args()
decals = DECALS(ns.conf)

confcut = {
    'z': ns.sigma_z,
    'r': ns.sigma_r,
    'g': ns.sigma_g,
    }

def apply_confcut(confidence, confcut):
    result = np.ones(len(confidence), dtype='?')
    for band in confcut:
        iband = product.bands[band]
        result &= confidence[:, iband] > confcut[band]
    return result


with h5py.File(ns.catalogue, 'r') as ff:
    mask = np.ones(len(ff['RA'][:]), dtype='?')

    if 'CONFIDENCE' in ff:
        confmask = apply_confcut(ff['CONFIDENCE'][:], confcut)
        mask &= confmask

    if ns.use_tycho_veto is not None:
        vetomask = ~ff['TYCHO_VETO'][ns.use_tycho_veto]
        mask &= vetomask

    ra = ff['RA'][:][mask]
    dec = ff['DEC'][:][mask]
    fc = ff['COMPLETENESS'][:][mask]
    h = ['RA', 'DEC', 'COMPETENESS']
    l = [ra, dec, fc]

    for band in ns.bands:
        iband = product.bands[band]
        if 'INTRINSIC_FLUX' in ff:
            h.append(band + 'flux')
            l.append(ff['INTRINSIC_FLUX'][:, iband][mask])
        h.append(band + 'noise')
        l.append(ff['INTRINSIC_NOISELEVEL'][:, iband][mask])

with open(ns.output, 'wb') as ff:
    ff.write(bytes('\t'.join(h) + '\n', encoding='utf8'))
    np.savetxt(ff, np.array(l).T, fmt='%12f')
