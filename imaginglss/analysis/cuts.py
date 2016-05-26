# "Helper" routines to centralize the flux and color cuts
# in one place where they can be called from other routines.

import numpy as N

__author__ = "Yu Feng and Martin White"
__version__ = "1.0"
__email__  = "yfeng1@berkeley.edu or mjwhite@lbl.gov"

def apply(comm, query, data):
    """ Apply a query, or 'cut' to data. 

        The Queres
    """
    total = sum(comm.allgather(len(data)) )
    if comm.rank == 0: 
        print(total)

    for expr in query:
        if comm.rank == 0:
            print("Applying Cut :", expr)
        comm.barrier()
        mask = expr.apply(data)
        selected = sum(comm.allgather(mask.sum())) 
        if comm.rank == 0:
            print("%s : %d / %d = %g"
                % (
                str(expr), selected, total,
                1. * selected / total))
    mask = query.apply(data)
    selected = sum(comm.allgather(mask.sum())) 
    if comm.rank == 0:
        print("%s : %d / %d = %g"
            % (
            str(query), selected, total,
            1. * selected / total))
    return mask

