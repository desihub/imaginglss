import numpy
from decals import DECALS

from matplotlib.figure import Figure
from matplotlib.backends.backend_agg import FigureCanvasAgg
from matplotlib.gridspec import GridSpec
#mag = 22.5 - 2.5 * log10(getdepth(398599, 'z')[0].data)
def lsst_iband_model(i):
    """ equation 3.7 in LSST Science book (SB_3.pdf)
        cummulative number of galaxies as i band mag, per 
        arcmin**2
    """
    Ngal = 46. * 10. **(0.31 * (i - 25.))
    return Ngal

if __name__ == '__main__':
    decals = DECALS('.')

    print 'total number of objects', decals.galaxies.shape
    print 'current survey area', decals.observed_area
    dec = decals.galaxies['DEC']
    ra = decals.galaxies['RA']

    g = 22.5 - 2.5 * numpy.log10(decals.galaxies['DECAM_FLUX'][:, decals.bands['g']])
    r = 22.5 - 2.5 * numpy.log10(decals.galaxies['DECAM_FLUX'][:, decals.bands['r']])
    z = 22.5 - 2.5 * numpy.log10(decals.galaxies['DECAM_FLUX'][:, decals.bands['z']])
    
    #numpy.savetxt('flux.txt', zip(ra, dec, g, r, z), header="ra dec g r z")

    random = numpy.loadtxt('random.txt', dtype=[
            ('RA', 'f4'),
            ('DEC', 'f4'),
            ('dg', 'f4'),
            ('dr', 'f4'),
            ('dz', 'f4'),
            ])

    random['dg'] = 22.5 - 2.5 * numpy.log10(random['dg'])
    random['dr'] = 22.5 - 2.5 * numpy.log10(random['dr'])
    random['dz'] = 22.5 - 2.5 * numpy.log10(random['dz'])

    fig = Figure()
    
    survey_area = decals.observed_area * 3600
    
    def histband(m, random, survey_area, ax, ax2, label):
        irange = numpy.linspace(14, 22, 10)
        h, bins = numpy.histogram(m, range=(14, 22))

        h = h.cumsum()
        eff_area = survey_area * numpy.array([
                1.0 * (random > bins[i]).sum() / len(random)
                for i in range(1, len(bins))])
        
        ax.plot(bins[1:], h / eff_area)
        ax.plot(irange, lsst_iband_model(irange))
        ax.set_xlabel(label)
        ax.set_ylabel('Number Density per arcmin**2')
        ax.set_title('type != S')
        ax.set_yscale('log')
        ax2.plot(bins[1:], eff_area/ survey_area)
        ax2.set_ylabel('Covered ratio')
        ax2.set_ylim(0, 1.1)
    gs = GridSpec(2, 3)
    ax = fig.add_subplot(gs[0, 0])
    ax1 = fig.add_subplot(gs[1, 0])
    histband(g, random['dg'], survey_area, ax, ax1, 'g')
    ax = fig.add_subplot(gs[0, 1])
    ax1 = fig.add_subplot(gs[1, 1])
    histband(r, random['dr'], survey_area, ax, ax1, 'r')
    ax = fig.add_subplot(gs[0, 2])
    ax1 = fig.add_subplot(gs[1, 2])
    histband(z, random['dz'], survey_area, ax, ax1, 'z')

    canvas = FigureCanvasAgg(fig)
    fig.savefig("grz-hist.png")

    fig = Figure()
    mask = random['dg'] > 19.5
    ax = fig.add_subplot(111)
    ax.hist2d(random['RA'][mask], random['DEC'][mask])
    ax.set_xlabel('RA')
    ax.set_ylabel('DEC')
    ax.set_title('dg > 19.5')
    canvas = FigureCanvasAgg(fig)
    fig.savefig("g-mask-hist.png")
    
    fig = Figure()
    mask = random['dr'] > 18.5
    ax = fig.add_subplot(111)
    ax.hist2d(random['RA'][mask], random['DEC'][mask])
    ax.set_xlabel('RA')
    ax.set_ylabel('DEC')
    ax.set_title('dr > 18.5')
    canvas = FigureCanvasAgg(fig)
    fig.savefig("r-mask-hist.png")

    fig = Figure()
    ax = fig.add_subplot(111)
    ax.hist2d(ra, dec)
    ax.set_xlabel('RA')
    ax.set_ylabel('DEC')
    ax.set_title('type != S')
    canvas = FigureCanvasAgg(fig)
    fig.savefig("ra-dec-hist.png")
    
