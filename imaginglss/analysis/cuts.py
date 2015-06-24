# "Helper" routines to centralize the flux and color cuts
# in one place where they can be called from other routines.
# Contains two useful classes:
#	Fluxes --	contains methods which return True when
#			an object's fluxes pass a certain cut.
#	Completeness --	contains methods which return True when	a given
#			depth puts it in the 100% complete sample.
# as well as a routine to compute N-sigma flux limits (findlim).

import numpy as N

__author__ = "Yu Feng and Martin White"
__version__ = "1.0"
__email__  = "yfeng1@berkeley.edu or mjwhite@lbl.gov"

def apply(comm, query, data):
    total = sum(comm.allgather(len(data)) )
    for expr in query:
        mask = expr.visit(data)
        selected = sum(comm.allgather(mask.sum())) 
        if comm.rank == 0:
            print("%s : %d / %d = %g"
                % (
                str(expr), selected, total,
                1. * selected / total))
    mask = query.visit(data)
    selected = sum(comm.allgather(mask.sum())) 
    if comm.rank == 0:
        print("%s : %d / %d = %g"
            % (
            str(query), selected, total,
            1. * selected / total))
    return mask

