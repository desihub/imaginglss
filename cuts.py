
class Color:
    @staticmethod
    def LRG(rflux, zflux, wflux):
        raise UnimplementedError

    @staticmethod
    def ELG(rflux, gflux, zflux):
        mask = rflux > 10**((22.5-23.4) / 2.5)
        mask &= zflux > 10**((0.3) / 2.5) * rflux
        mask &= zflux < 10**((1.5) / 2.5) * rflux
        mask &= rflux ** 2< gflux * zflux * 10 ** (-0.2/2.5)
        mask &= zflux > gflux * 10**(1.2/2.5)
        return mask
 
class Completeness:
    """ Completeness cuts in intrinsic limiting flux 
        (nominal limit)
       
        for any band depth and MW transmission:

        lim = sigma / (dep) ** 0.5 / tran
    """
    @staticmethod
    def LRG(rlim, zlim, wlim):
        mask  = (rlim<10.0**((22.5-23.00)/2.5))
        mask &= (zlim<10.0**((22.5-20.56)/2.5))
        mask &= (wlim<10.0**((22.5-19.50)/2.5))
        return mask

    @staticmethod
    def ELG(rlim, glim, zlim):
        mask  = rlim < 10 ** ((22.5 - 23.4) / 2.5) 
        mask &= zlim < 10 ** ((22.5 - 23.4 + 0.3) / 2.5)
        mask &= glim < 10 ** ((22.5 - 23.4 - 1.5 + 0.2) / 2.5)
        return mask

    @staticmethod
    def QSO(rlim):
        mask  = rlim<10.0**((22.5-23.00)/2.5)
        return mask 

