from model.brickindex import BrickIndex
from model.datarelease import DataRelease

dr = DataRelease(root='.', version='EDR')
bi = dr.brickindex
ob = dr.observed_bricks

brick = bi[ob[0]]
print dr.images['DEPTH'].open(brick, band='z')
print brick
