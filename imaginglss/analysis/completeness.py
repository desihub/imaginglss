from __future__ import print_function
import numpy as np
from imaginglss.analysis    import targetselection

def CompletenessEstimator(datarelease, objecttype,
    fluxes, noises, confidence={'r': 5.0, 'g': 5.0, 'z': 3.0}):
    """
        Create a completeness estimator for an object type,
        based on the object type, intrinsic fluxes and intrinsic noise levels.

        Parameters
        ----------
        datarelease : 
            a data release object, where we look up the translation between
            band name and band index.
        objectype : string
            type of object. We use this string to look up the target definition
            and the involved bands in the definition.
        fluxes : array_like (N, 6)
            intrinsic fluxes from DECAM; usually calculated by imglss-mpi-select-objects.py.
            in nano-maggies.
        noises : array_like (N, 6)
            intrinsic 1-sigma noise level from DECAM; usually calculated by imglss-mpi-select-objects.py,
            in nano-maggies.
        confidence : dictionary, optional
            confidence (in sigma) for each band.
    """

    # we use the integrator in kdcount for the estimator
    from kdcount import KDTree

    # the targetselection type knows which bands are useful
    # we only put limits on those bands.
    #
    # FIXME: @ekitanidis what about adding a link to your talk slides 
    # explaining this?

    fluxcut = getattr(targetselection, objecttype)

    fluxes = np.array([
        fluxes[:, datarelease.bands[band]]
        for band in fluxcut.bands]).T

    noises = np.array([
        confidence[band] * noises[:, datarelease.bands[band]]
        for band in fluxcut.bands]).T

    # This will be the 100% completeness limit for the given confidence
    lim = fluxes.min(axis=0)

    mask = (noises <= lim).all(axis=-1)
    model = fluxes[mask]
    tree = KDTree(model)
    root = tree.root

    def fcmodelfunc(noises):
        noises = np.array([
            confidence[band] * noises[:, datarelease.bands[band]]
            for band in fluxcut.bands]).T
        seen = root.integrate(noises, np.inf)
        mask = (noises <= lim).all(axis=-1)

        # Watchout:
        # Only 100% complete area has fcomp == 1.0
        # otherwise we give a completeness slightly less than 1.0

        fcomp = 1.0 * seen / (len(model) + 1.0)
        fcomp[mask] = 1.0
        return fcomp

    return fcmodelfunc
