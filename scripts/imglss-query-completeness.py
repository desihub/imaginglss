
from __future__ import print_function
import numpy as np
import h5py

from imaginglss             import DECALS
from imaginglss.model             import dataproduct
from imaginglss.analysis    import targetselection
from imaginglss.analysis    import completeness
from imaginglss.utils       import output
from imaginglss.analysis    import tycho_veto

from argparse import ArgumentParser

ap = ArgumentParser()
ap.add_argument("ObjectType", choices=[i for i in targetselection.__all__])

ap.add_argument("objects", 
        help="object catalogue for building the completeness model.")
ap.add_argument("query", 
        help="catalogue to query completeness")

allvetos = [i for i in dir(tycho_veto) if not str(i).startswith( '_' )]
ap.add_argument("--use-tycho-veto", type=str, choices=allvetos, default=None)
ap.add_argument("--sigma-z", type=float, default=3.0)
ap.add_argument("--sigma-g", type=float, default=5.0)
ap.add_argument("--sigma-r", type=float, default=5.0)

ap.add_argument("--conf", default=None,
        help="Path to the imaginglss config file, default is from DECALS_PY_CONFIG")

ns = ap.parse_args()
decals = DECALS(ns.conf)

np.seterr(divide='ignore', invalid='ignore')

def query_completeness(decals, ns):
    with h5py.File(ns.objects, 'r') as ff:
        FLUX = ff['INTRINSIC_FLUX'][:]
        NOISE = ff['INTRINSIC_NOISELEVEL'][:]
        if ns.use_tycho_veto is not None:
            veto = ff['TYCHO_VETO'][ns.use_tycho_veto][:]
            FLUX = FLUX[~veto]
            NOISE = NOISE[~veto]

    confidence = {
        'z': ns.sigma_z,
        'r': ns.sigma_r,
        'g': ns.sigma_g,
        }

    bands = dataproduct.bands
    # a list of limit_bands in integer indices
    active_bands = getattr(targetselection, ns.ObjectType).limit_bands

    active_ibands = [bands[band] for band in active_bands]

    active_conf  = np.array([ confidence[band] for band in active_bands])

    # if a band is used and unvisited we will exclude it from
    # the model, and later require the completeness as 0.
    # (it is outside of survey for this target type)

    all_bands = getattr(targetselection, ns.ObjectType).bands
    all_ibands = np.array([ bands[band] for band in all_bands])

    visitedmask = (~np.isinf(NOISE[:, all_ibands])).all(axis=-1)
    print ('Number of visited objects', visitedmask.sum())
    active_flux = FLUX[:, active_ibands][visitedmask]

    active_noise = NOISE[:, active_ibands][visitedmask]

    model = completeness.CompletenessEstimator(active_flux, active_noise, active_conf)

    with h5py.File(ns.query, 'r') as ff:
        NOISE = ff['INTRINSIC_NOISELEVEL'][:]

    FC = model(NOISE[:, active_ibands])
    visitedmask = (~np.isinf(NOISE[:, all_ibands])).all(axis=-1)
    FC[~visitedmask] = 0

    return FC

if __name__=="__main__":

    FC = query_completeness(decals, ns)
    print("N = %d Max FC = %g, MIN FC= %g" % (len(FC), FC.max(), FC.min()))
    with h5py.File(ns.query, 'r+') as ff:
        if 'COMPLETENESS' in ff:
            del ff['COMPLETENESS']
        ds = ff.create_dataset('COMPLETENESS', data=FC)
        for k, v in ns.__dict__.items():
            ds.attrs[k] = v
