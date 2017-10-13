import h5py
import numpy as np

from   imaginglss             import DECALS
from imaginglss.model import dataproduct

from imaginglss.cli import CLI

cli = CLI("Naive Angular Cross Correlation Function",
        enable_target_plugins=True,
        enable_tycho_veto=True,
        enable_confidence=True)
cli.add_argument("--np", default=8, type=int, help="nubmer of in-node cores to use")
cli.add_argument("catalogue1", help="internal catalogue of HDF5 type.")
cli.add_argument("random1", help="internal catalogue of HDF5 type.")
cli.add_argument("catalogue2", help="internal catalogue of HDF5 type.")
cli.add_argument("output", help="text file to store the correlation function." )

from kdcount import sphere
from kdcount import correlate

ns = cli.parse_args()

def corr():
    data1file = h5py.File(ns.catalogue1, 'r')
    data2file = h5py.File(ns.catalogue2, 'r')
    rand1file = h5py.File(ns.random1, 'r')

    data1mask = data1file['COMPLETENESS'][:] >= 1
    for vetoname in ns.use_tycho_veto:
        data1mask &= ~data1file['TYCHO_VETO'][vetoname][:]
    data1RA = data1file['RA'][:][data1mask]
    data1DEC = data1file['DEC'][:][data1mask]

    rand1mask = rand1file['COMPLETENESS'][:] >= 1
    for vetoname in ns.use_tycho_veto:
        rand1mask &= ~rand1file['TYCHO_VETO'][vetoname][:]
    rand1RA = rand1file['RA'][:][rand1mask]
    rand1DEC = rand1file['DEC'][:][rand1mask]

#    data2mask = data2file['COMPLETENESS'][:] >= 1
#    for vetoname in ns.use_tycho_veto:
#        data2mask &= ~data2file['TYCHO_VETO'][vetoname][:]
    data2mask = Ellipsis
    data2RA = data2file['RA'][:][data2mask]
    data2DEC = data2file['DEC'][:][data2mask]

    data1 = sphere.points(data1RA, data1DEC)
    data2 = sphere.points(data2RA, data2DEC)
    rand1 = sphere.points(rand1RA, rand1DEC)
    abin = sphere.AngularBinning(np.logspace(-3, 0, 16, endpoint=True))

    DD = correlate.paircount(data1, data2, abin, np=ns.np)
    DR = correlate.paircount(rand1, data2, abin, np=ns.np)

    r =  1. * len(data1) / len(rand1)

    dd = 1.0 * DD.sum1
    dr = 1.0 * DR.sum1 * r

    return abin.angular_centers, (dd - dr)/ dr


if __name__ == '__main__':

    t, w = corr()
    np.savetxt(ns.output, np.array([t, w]).T, header="theta w")

