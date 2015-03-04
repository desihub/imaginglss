from model.datarelease import DataRelease

dr = DataRelease()

flux = dr.catalogue['DECAM_FLUX'].T
RA = dr.catalogue['RA']
DEC = dr.catalogue['DEC']

print flux
primary = dr.catalogue['BRICK_PRIMARY']
print primary
extinction = dr.catalogue['EXTINCTION'].T
print extinction

GFLUX = flux[1] * 10 ** (extinction[1] / 2.5)
RFLUX = flux[2] * 10 ** (extinction[2] / 2.5)
ZFLUX = flux[4] * 10 ** (extinction[4] / 2.5)
mask = primary == 1
print mask.shape, RFLUX.shape
mask &= RFLUX > 10**((22.5-23.4) / 2.5)
mask &= ZFLUX > 10**((0.3) / 2.5) * RFLUX
mask &= ZFLUX < 10**((1.5) / 2.5) * RFLUX
mask &= RFLUX ** 2< GFLUX * ZFLUX * 10 ** (-0.2/2.5)
mask &= ZFLUX > GFLUX * 10**(1.2/2.5)


print 'selected', mask.sum(), 'items'
