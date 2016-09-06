import h5py

class WISE(numpy.ndarray):
    """ Representing a WISE catalogue that
        is used to veto objects near stars 
    """
    def __new__(kls, path):
        data = h5py.File(path, 'r')
        self = numpy.empty(len(data['RA']),
            dtype=[
                ('RA', 'f8'),
                ('DEC', 'f8'),
                ('W1MPRO', 'f8'),
                ])
        self['RA'] = data['RA']
        self['DEC'] = data['DEC']
        self['W1MPRO'] = data['W1MPRO']
        return self.view(type=kls)
