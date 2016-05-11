from __future__ import print_function
import numpy as np
from imaginglss.analysis    import targetselection

def CompletenessEstimator(fluxes, noises, confidence):
    """
        Create a completeness estimator for an object type,
        based on the object type, intrinsic fluxes and intrinsic noise levels.

        Parameters
        ----------
        fluxes : array_like (N, Nbands)
            intrinsic fluxes from DECAM; usually calculated by imglss-mpi-select-objects.py.
            in nano-maggies.
        noises : array_like (N, Nbands)
            intrinsic 1-sigma noise level from DECAM; usually calculated by imglss-mpi-select-objects.py,
            in nano-maggies.
        confidence : array_like (Nbands)
            confidence (in sigma) for each band.
    """

    # we use the integrator in kdcount for the estimator
    from kdcount import KDTree

    # the targetselection type knows which bands are useful
    # we only put limits on those bands.
    #
    # FIXME: @ekitanidis what about adding a link to your talk slides 
    # explaining this?

    # This will be the 100% completeness limit for the given confidence
    lim = fluxes.min(axis=0)

    noises = confidence[None, :] * noises
    mask = (noises <= lim).all(axis=-1)
    model = fluxes[mask]
    tree = KDTree(model)
    root = tree.root

    def fcmodelfunc(query_noises):
        query_noises = confidence[None, :] * query_noises
        seen = root.integrate(query_noises, np.inf)
        mask = (query_noises <= lim).all(axis=-1)

        # Watchout:
        # Only 100% complete area has fcomp == 1.0
        # otherwise we give a completeness slightly less than 1.0

        fcomp = 1.0 * seen / (len(model) + 1.0)
        fcomp[mask] = 1.0
        return fcomp

    return fcmodelfunc
