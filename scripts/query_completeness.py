
from __future__ import print_function
import numpy as np

from imaginglss             import DECALS
from imaginglss.analysis    import targetselection
from imaginglss.utils       import output
from kdcount import KDTree

from argparse import ArgumentParser

ap = ArgumentParser("query_completeness.py")
ap.add_argument("ObjectType", choices=[i for i in targetselection.__all__])
ap.add_argument("output", type=output.writer)
ap.add_argument("objects", type=output.writer)
ap.add_argument("noises", type=output.writer)
ap.add_argument("--conf", default=None,
        help="Path to the imaginglss config file, default is from DECALS_PY_CONFIG")

ns = ap.parse_args()
ns.conf = DECALS(ns.conf)

np.seterr(divide='ignore', invalid='ignore')

def build_model(ns, fluxes, noises, sigmas={'r':3, 'g':3, 'z':3}):
    fluxcut = getattr(targetselection, ns.ObjectType)
    
    fluxes = np.array([
        fluxes['DECAM_INTRINSIC_FLUX'][:, ns.conf.datarelease.bands[band]]
        for band in fluxcut.bands]).T

    noises = np.array([
        sigmas[band] * noises['DECAM_INTRINSIC_NOISE_LEVEL'][:, ns.conf.datarelease.bands[band]]
        for band in fluxcut.bands]).T

    mask = (fluxes > noises).all(axis=-1)
    model = fluxes[mask]
    tree = KDTree(model)
    root = tree.root 

    def modelfunc(noises):
        noises = np.array([
            sigmas[band] * noises['DECAM_INTRINSIC_NOISE_LEVEL'][:, ns.conf.datarelease.bands[band]]
            for band in fluxcut.bands]).T
        seen = root.integrate(noises, np.inf)
        
        return 1.0 * seen / len(model)

    return modelfunc

def query_completeness(ns):
    object_fluxes = ns.objects.read('FLUXES')
    object_noises = ns.objects.read('NOISES')
    model = build_model(ns, object_fluxes, object_noises)

    noises = ns.noises.read('NOISES')

    dtype = [('FRACTION_COMPLETENESS', 'f8')]
    FC = np.empty(len(noises), dtype=dtype)
    FC['FRACTION_COMPLETENESS'][:] = model(noises) 
    return FC

if __name__=="__main__":

    FC = query_completeness(ns)

    ns.output.write(FC, ns.__dict__, 'FC')
