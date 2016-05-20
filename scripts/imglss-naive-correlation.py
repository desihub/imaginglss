import h5py
import numpy as np

from   imaginglss             import DECALS
from imaginglss.model import dataproduct

from imaginglss.cli import CLI

cli = CLI("Naive Angular Correlation Function",
        enable_target_plugins=True,
        enable_tycho_veto=True,
        enable_confidence=True)
cli.add_argument("--np", default=8, type=int, help="nubmer of in-node cores to use")
cli.add_argument("catalogue", help="internal catalogue of HDF5 type.")
cli.add_argument("random", help="internal catalogue of HDF5 type.")
cli.add_argument("output", help="text file to store the correlation function." )

from kdcount import sphere
from kdcount import correlate

ns = cli.parse_args()

def corr():
    datafile = h5py.File(ns.catalogue, 'r')
    randfile = h5py.File(ns.random, 'r')

    datamask = datafile['COMPLETENESS'][:] >= 1
    if ns.use_tycho_veto is not None:
        datamask &= ~datafile['TYCHO_VETO'][ns.use_tycho_veto][:]
    dataRA = datafile['RA'][:][datamask]
    dataDEC = datafile['DEC'][:][datamask]

    randmask = randfile['COMPLETENESS'][:] >= 1
    if ns.use_tycho_veto is not None:
        randmask &= ~randfile['TYCHO_VETO'][ns.use_tycho_veto][:]
    randRA = randfile['RA'][:][randmask]
    randDEC = randfile['DEC'][:][randmask]

    data = sphere.points(dataRA, dataDEC)
    rand = sphere.points(randRA, randDEC)
    abin = sphere.AngularBinning(np.logspace(-3, 0, 16, endpoint=True))

    DD = correlate.paircount(data, data, abin, np=ns.np)
    DR = correlate.paircount(data, rand, abin, np=ns.np)
    RR = correlate.paircount(rand, rand, abin, np=ns.np)

    r =  1. * len(data) / len(rand)

    dd = 1.0 * DD.sum1
    dr = 1.0 * DR.sum1 * r
    rr =  1.0 * RR.sum1 * (r * r)

    return abin.centers, (dd - 2 * dr + rr) / rr


if __name__ == '__main__':

    t, w = corr()
    np.savetxt(ns.output, zip(t, w), header="theta w")

