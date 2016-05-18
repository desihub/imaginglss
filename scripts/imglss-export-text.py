import h5py
import numpy as np

from   imaginglss             import DECALS
from imaginglss.model import dataproduct

from imaginglss.cli import CLI

cli = CLI("Export Final Data Products",
        enable_target_plugins=True,
        enable_tycho_veto=True,
        enable_confidence=True)

cli.add_target_type_argument("ObjectType")
cli.add_argument("catalogue", help="internal catalogue of HDF5 type.")
cli.add_argument("output", help="text file to store the catalogue." )
cli.add_argument("--bands", nargs='+', type=str, 
        choices=dataproduct.bands.keys(), 
        default=['r', 'g', 'z'], help="Bands to export in the text file")

ns = cli.parse_args()
decals = DECALS(ns.conf)

confcut = {
    'z': ns.sigma_z,
    'r': ns.sigma_r,
    'g': ns.sigma_g,
    }

def apply_confcut(objecttype, confidence, confcut):
    result = np.ones(len(confidence), dtype='?')
    for band in confcut:
        if not band in objecttype.mag_bands: continue
        iband = dataproduct.bands[band]
        result &= confidence[:, iband] > confcut[band]
    return result


with h5py.File(ns.catalogue, 'r') as ff:
    mask = np.ones(len(ff['RA'][:]), dtype='?')

    if 'CONFIDENCE' in ff:
        # this is a catalogue, apply confidence cut.
        confmask = apply_confcut(ns.ObjectType, ff['CONFIDENCE'][:], confcut)
        mask &= confmask

    if ns.use_tycho_veto is not None:
        if not 'TYCHO_VETO' in ff:
            raise KeyError("TYCHO_VETO dataset is not found. Run imglss-query-tycho-veto on `%s` first." % ns.catalogue)
        vetomask = ~ff['TYCHO_VETO'][ns.use_tycho_veto]
        mask &= vetomask

    ra = ff['RA'][:][mask]
    dec = ff['DEC'][:][mask]

    if not 'COMPLETENESS' in ff:
        raise KeyError("COMPLETENESS dataset is not found. Run imglss-query-completeness on `%s` first." % ns.catalogue)

    fc = ff['COMPLETENESS'][:][mask]
    h = ['RA', 'DEC', 'COMPETENESS']
    l = [ra, dec, fc]

    for band in ns.bands:
        iband = dataproduct.bands[band]
        if 'INTRINSIC_FLUX' in ff:
            # this is a catalogue, dump flux as well as noise
            h.append(band + 'flux')
            l.append(ff['INTRINSIC_FLUX'][:, iband][mask])
        h.append(band + 'noise')
        l.append(ff['INTRINSIC_NOISELEVEL'][:, iband][mask])

with open(ns.output, 'wb') as ff:
    ff.write(('# ' + '\t'.join(h) + '\n').encode())
    np.savetxt(ff, np.array(l).T, fmt='%12f')
