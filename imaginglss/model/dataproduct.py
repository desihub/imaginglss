import numpy

bands = {'u':0, 'g':1, 'r':2, 'i':3, 'z':4, 'Y':5, 'W1': 6, 'W2' : 7, 'W3' : 8, 'W4' : 9}

Nbands = len(bands)

ObjectCatalogue = numpy.dtype(
    [
        ('RA', 'f8'),
        ('DEC', 'f8'),
        ('INTRINSIC_FLUX', ('f8', Nbands)),
        ('INTRINSIC_NOISELEVEL', ('f8', Nbands)),
        ('CONFIDENCE', ('f8', Nbands)),
    ])

RandomCatalogue = numpy.dtype(
    [
        ('RA', 'f8'),
        ('DEC', 'f8'),
        ('INTRINSIC_NOISELEVEL', ('f8', Nbands)),
    ])

