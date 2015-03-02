from model.brickindex import BrickIndex
from model.datarelease import DataRelease
import numpy

dr = DataRelease(root='.', version='EDR')
bi = dr.brickindex
ob = dr.observed_bricks

brick = bi[ob[0]]
print dr.images['DEPTH'].open(brick, band='z')
print brick

def test398599():
    """ test image readout on brick-398599. 
    """
    dr = DataRelease(root='.', version='EDR')
    bi = dr.brickindex
    ob = dr.observed_bricks

    brick = bi[ob[0]]

    print ob[0], brick
    x, y = numpy.indices((3600, 3600))
    x = numpy.ravel(x) + 0.5
    y = numpy.ravel(y) + 0.5
    coord = brick.revert((x, y))
    img = dr.readout(coord, keys=(('IMAGE', 'z'), ))[..., 0]

    print 'found', (~numpy.isnan(img)).sum()
    img2 = dr.images['IMAGE'].open(brick, band='z')
    diff = img.reshape(3600, 3600) - img2

    # FIXME: tighten this up
    assert (diff[300:-300, 300:-300] == 0).all()
    print 'passed'

def testcat():
    dr = DataRelease(root='.', version='EDR')
    print dr.catalogue 

testcat()
test398599()
