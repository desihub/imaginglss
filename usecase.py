from model.brickindex import BrickIndex
from model.datarelease import DataRelease
import numpy

dr = DataRelease(version='EDR')
bi = dr.brickindex
ob = dr.observed_brickids

brick = bi[ob[0]]
print dr.images['z']['IMAGE'].metadata(brick)
#print dr.images['DEPTH'].open(brick, band='z')
#print brick

def test398599():
    """ test image readout on brick-398599. 
    """
    dr = DataRelease(version='EDR')
    bi = dr.brickindex
    ob = dr.observed_brickids
    brick = bi[ob[0]]
    dr.images['z']['IMAGE'].preload([brick])

    img2 = dr.images['z']['IMAGE'].open(brick)
    assert (img2[300:-300, 300:-300] != 0).any()

    print ob[0], brick
    x, y = numpy.indices((3600, 3600))
    x = numpy.ravel(x) + 0.5
    y = numpy.ravel(y) + 0.5
    coord = brick.revert((x, y))

    img3 = brick.readout(coord, dr.images['z']['IMAGE'])
    diff = img3.reshape(3600, 3600) - img2
    assert (diff[300:-300, 300:-300] == 0).all()

    x2, y2 = brick.query(coord)
    assert numpy.allclose(x, x2)
    assert numpy.allclose(y, y2)


    img = dr.readout(coord, (dr.images['z']['IMAGE'],))[..., 0]
    print 'brick readout passed'
    print 'found', (~numpy.isnan(img)).sum()
    diff = img.reshape(3600, 3600) - img2
    assert (diff[300:-300, 300:-300] == 0).all()
    print 'dr readout passed'


def testrandom():
    u1, u2 = numpy.random.random(size=(2, 4000000))
    RA = (246 - 239) * u1 + 239
    a = 0.5 * ((numpy.cos(12. / 180. * numpy.pi)  + 1))
    b = 0.5 * ((numpy.cos(5. / 180. * numpy.pi)  + 1))
    u2 = (a - b) * u2 + b
    
    DEC = numpy.arccos(2.0 * u2 - 1) * (180. / numpy.pi)

    dr = DataRelease(version='EDR')

    coord = (RA, DEC)
    depth = dr.readout(coord, (dr.images['u']['IMAGE'],))[..., 0]
    print depth
 
def testcat():
    dr = DataRelease(version='EDR')
    print dr.catalogue 
    print dr.catalogue['RA']
#testrandom()
testcat()
test398599()
