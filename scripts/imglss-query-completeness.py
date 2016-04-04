
from __future__ import print_function
import numpy as np

from imaginglss             import DECALS
from imaginglss.analysis    import targetselection
from imaginglss.analysis    import completeness
from imaginglss.utils       import output

from argparse import ArgumentParser

ap = ArgumentParser()
ap.add_argument("ObjectType", choices=[i for i in targetselection.__all__])

ap.add_argument("output", type=output.writer,
        help="Output file name. will write to the FC extension")
ap.add_argument("objects", type=output.writer,
        help="Object file name to build the model. NOISES and FLUXES extensions are used")
ap.add_argument("noises", type=output.writer,
        help="Input file name to query the completeness for NOISES extension is used.")

ap.add_argument("--sigma-z", type=float, default=3.0)
ap.add_argument("--sigma-g", type=float, default=5.0)
ap.add_argument("--sigma-r", type=float, default=5.0)

ap.add_argument("--conf", default=None,
        help="Path to the imaginglss config file, default is from DECALS_PY_CONFIG")

ns = ap.parse_args()
ns.conf = DECALS(ns.conf)

np.seterr(divide='ignore', invalid='ignore')


def query_completeness(ns):
    object_fluxes = ns.objects.read('FLUXES')['DECAM_INTRINSIC_FLUX']
    object_noises = ns.objects.read('NOISES')['DECAM_INTRINSIC_NOISE_LEVEL']
    confidence = {
        'z': ns.sigma_z,
        'r': ns.sigma_r,
        'g': ns.sigma_g,
        }
    model = completeness.CompletenessEstimator(ns.conf.datarelease, ns.ObjectType, object_fluxes, object_noises, confidence)

    noises = ns.noises.read('NOISES')['DECAM_INTRINSIC_NOISE_LEVEL']

    dtype = [ ('FRACTION_COMPLETENESS', 'f8')]
    FC = np.empty(len(noises), dtype=dtype)
    FC['FRACTION_COMPLETENESS'][:] = model(noises)
    return FC

if __name__=="__main__":

    FC = query_completeness(ns)

    ns.output.write(FC, ns.__dict__, 'FC')
