from model.datarelease import DataRelease
import numpy

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
mask &= RFLUX > 10**((22.5-23.4) / 2.5)
mask &= ZFLUX > 10**((0.3) / 2.5) * RFLUX
mask &= ZFLUX < 10**((1.5) / 2.5) * RFLUX
mask &= RFLUX ** 2< GFLUX * ZFLUX * 10 ** (-0.2/2.5)
mask &= ZFLUX > GFLUX * 10**(1.2/2.5)

print 'selected', mask.sum(), 'items', 'out of', len(mask)
go = 22.5 - 2.5 * numpy.log10(flux[1][mask])
ro = 22.5 - 2.5 * numpy.log10(flux[2][mask])
zo = 22.5 - 2.5 * numpy.log10(flux[4][mask])
GFLUX = flux[1] * 10 ** (extinction[1] / 2.5)
RFLUX = flux[2] * 10 ** (extinction[2] / 2.5)
ZFLUX = flux[4] * 10 ** (extinction[4] / 2.5)

gi = 22.5 - 2.5 * numpy.log10(GFLUX[mask])
ri = 22.5 - 2.5 * numpy.log10(RFLUX[mask])
zi = 22.5 - 2.5 * numpy.log10(ZFLUX[mask])

output = numpy.empty(mask.sum(), dtype=[
        ('RA', 'f8'),
        ('DEC', 'f8'),
        ('gi', 'f4'),
        ('ri', 'f4'),
        ('zi', 'f4'),
        ('go', 'f4'),
        ('ro', 'f4'),
        ('zo', 'f4'),
        ])
output['RA'] = RA[mask]
output['DEC'] = DEC[mask]
output['gi'] = gi
output['ri'] = ri
output['zi'] = zi
output['go'] = go
output['ro'] = ro
output['zo'] = zo
numpy.savetxt('elg-targets.txt', output, header='#RA, DEC gi ri zi go ro zo')

from matplotlib.figure import Figure
from matplotlib.backends.backend_agg import FigureCanvasAgg

fig = Figure()
ax = fig.add_subplot(111)
ax.hist2d(output['RA'], output['DEC'])
canvas = FigureCanvasAgg(fig)
fig.savefig('elg-ra-dec.png')

fig = Figure()
ax = fig.add_subplot(111)
ax.plot(output['ri'], output['ri'] - output['zi'], '. ')
ax.axvline(23.4, label='r < 0.3')
ax.axhline(0.3, label='(r-z) > 0.3')
ax.axhline(1.5, label='(r-z) < 1.5')
ax.set_xlabel('Intrinsic z')
ax.set_ylabel('Intrinsic r-z')
ax.grid()
canvas = FigureCanvasAgg(fig)
fig.savefig('elg-ri-rz.png')

fig = Figure()
ax = fig.add_subplot(111)
x = numpy.linspace(
        (output['gi'] - output['ri']).min(),
        (output['gi'] - output['ri']).max(), 100)

ax.plot(output['gi'] - output['ri'], output['ri'] - output['zi'], '. ')
ax.plot(x, x + 0.2, label='(g-r) < (r-z) - 0.2')
ax.plot(x, 1.2 - x , label='(g-r) > 1.2 - (r-z)')
ax.set_xlabel('Intrinsic g-r')
ax.set_ylabel('Intrinsic r-z')
ax.grid()
canvas = FigureCanvasAgg(fig)
fig.savefig('elg-gr-rz.png')
