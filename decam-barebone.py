import numpy
from decals import DECALS

from matplotlib.figure import Figure
from matplotlib.backends.backend_agg import FigureCanvasAgg

#mag = 22.5 - 2.5 * log10(getdepth(398599, 'z')[0].data)

if __name__ == '__main__':
    decals = DECALS('.')

    print 'total number of objects', decals.galaxies.shape
    dec = decals.galaxies['DEC']
    ra = decals.galaxies['RA']

    g = decals.galaxies['DECAM_FLUX'][:, decals.bands['g']]
    r = decals.galaxies['DECAM_FLUX'][:, decals.bands['r']]
    z = decals.galaxies['DECAM_FLUX'][:, decals.bands['z']]
    
    #numpy.savetxt('flux.txt', zip(ra, dec, g, r, z), header="ra dec g r z")

    fig = Figure()
    ax = fig.add_subplot(131)
    ax.hist(22.5 - 2.5 * numpy.log10(g), range=(14, 22), log=True)
    ax.set_xlabel("g")
    ax.set_title('type != S')
    
    ax = fig.add_subplot(132)
    ax.hist(22.5 - 2.5 * numpy.log10(r), range=(14, 22), log=True)
    ax.set_xlabel("r")
    ax.set_title('type != S')
    ax = fig.add_subplot(133)
    ax.hist(22.5 - 2.5 * numpy.log10(z), range=(14, 22), log=True)
    ax.set_xlabel("z")
    ax.set_title('type != S')
    canvas = FigureCanvasAgg(fig)
    fig.savefig("grz-hist.png")

    fig = Figure()
    ax = fig.add_subplot(111)
    ax.hist2d(ra, dec)
    ax.set_xlabel('RA')
    ax.set_ylabel('DEC')
    ax.set_title('type != S')
    canvas = FigureCanvasAgg(fig)
    fig.savefig("ra-dec-hist.png")
    
