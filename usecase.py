from model.datarelease import DataRelease
from model.imagerepo import ImageRepo
from pprint import pprint
import numpy

dr = DataRelease()

brick = dr.observed_bricks[0]
print dr.images['image']['z'].metadata(brick)
print '\n'.join(sorted([str(f) for f in dr.catalogue.dtype.fields]))

def test398599():
    """ test image readout on brick-398599. 
    """
    dr = DataRelease()
    brick = dr.observed_bricks[0]
    dr.images['image']['z'].preload([brick])

    img2 = dr.images['image']['z'].open(brick)
    assert (img2[300:-300, 300:-300] != 0).any()

    print 'Testing on', brick
    print dr.images['image']['z'].get_filename(brick)
    y, x = numpy.indices(img2.shape)
    x = numpy.ravel(x) + 0.5
    y = numpy.ravel(y) + 0.5
    coord = brick.revert(dr.images['image']['z'], (x, y))

    x2, y2 = brick.query(dr.images['image']['z'], coord)
    print x, x2
    print y, y2
    assert numpy.allclose(x, x2)
    assert numpy.allclose(y, y2)


    img3 = brick.readout(coord, dr.images['image']['z'])
    diff = img3.reshape(img2.shape) - img2
    assert (diff[300:-300, 300:-300] == 0).all()


    img = dr.readout(coord, dr.images['image']['z'])
    print 'brick readout passed'
    print 'found', (~numpy.isnan(img)).sum()
    diff = img.reshape(img2.shape) - img2
    assert (diff[300:-300, 300:-300] == 0).all()
    print 'dr readout passed'


def testrandom():
    dr = DataRelease()
    ebv = dr.images['ebv']
    print dr.observed_range
    u1, u2 = numpy.random.random(size=(2, 4))
    RA = (dr.observed_range.ramax - dr.observed_range.ramin) * u1 + dr.observed_range.ramin
    a = 0.5 * ((numpy.cos(dr.observed_range.decmax / 180. * numpy.pi)  + 1))
    b = 0.5 * ((numpy.cos(dr.observed_range.decmin / 180. * numpy.pi)  + 1))
    u2 = (a - b) * u2 + b
    
    DEC = numpy.arccos(2.0 * u2 - 1) * (180. / numpy.pi)

#    print 'prefetched'
    coord = (RA, DEC)
    print coord
    depth = dr.readout(coord, dr.images['depth']['z'])
    ebv = dr.readout(coord, ebv)
    print depth, ebv

def testebv():
    dr = DataRelease()
    RA = dr.catalogue['RA']
    DEC = dr.catalogue['DEC']
    EXT_R = dr.catalogue['DECAM_EXTINCTION'][:, 2]
    ebv = dr.images['ebv']
    coord = (RA, DEC)
    ebv = dr.readout(coord, ebv)
    good = ~numpy.isnan(ebv)
    print 'matched', good.sum(), 'out of', len(RA)
    ebv = ebv[good] * 2.165
    EXT_R = EXT_R[good]
    print (ebv-EXT_R).min()
    assert numpy.allclose(ebv, EXT_R, 1e-2)

def testcat():
    dr = DataRelease()
    print dr.catalogue 
    print dr.catalogue['RA']
testebv()
testrandom()
testcat()
test398599()
