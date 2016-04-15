import numpy

def radec2pos(ra, dec):
    """ converting ra dec to position on a unit sphere.
        ra, dec are in degrees.
    """
    pos = numpy.empty(len(ra), dtype=('f8', 3))
    ra = ra * (numpy.pi / 180)
    dec = dec * (numpy.pi / 180)
    pos[:, 2] = numpy.sin(dec)
    pos[:, 0] = numpy.cos(dec) * numpy.sin(ra)
    pos[:, 1] = numpy.cos(dec) * numpy.cos(ra)
    return pos

def veto(coord, center, R):
    """
        Returns a veto mask for coord. any coordinate within R of center
        is vet.

        Parameters
        ----------
        coord : (RA, DEC)
        center : (RA, DEC)
        R     : degrees

        Returns
        -------
        Vetomask : True for veto, False for keep.

    """
    from kdcount import KDTree

    pos = radec2pos(center[0], center[1])
    tree = KDTree(pos)

    R = 2 * numpy.sin(numpy.radians(R) * 0.5)

    pos = radec2pos(coord[0], coord[1])
    other = KDTree(pos)
    vetoflag = numpy.zeros(len(coord), dtype='?')

    Rmax = R.max()

    def process(r, i, j):
        # i is tycho, j is objects
        rcut = R[i]
        jcut = j[r < rcut]
        vetoflag[jcut] |= True

    tree.root.enum(other.root, Rmax, process)
    return vetoflag
