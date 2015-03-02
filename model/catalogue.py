import numpy

class Catalogue(object):
    def __init__(self, data):
        self.data = data

        self.DEC = data['DEC']
        self.RA = data['RA']
        self.TYPE = data['TYPE']
        self.BRICKID = data['BRICKID']

    def neighbours(self, dec, ra, sep):
        pass
