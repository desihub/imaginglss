import numpy

class ColumnStore(object):
    """ a cached column store 
        subclass shall implement
    """
    def __init__(self, dtype=None):
        self.cache = {}
        self.dtype = dtype

    def __setitem__(self, column, value):
        value = numpy.array(value)
        self.cache[column] = value

        d = dict(self.dtype.fields)
        subshape = value.shape[1:]
        if len(subshape) == 0:
            d[column] = (value.dtype.base, ())
        else:
            d[column] = (value.dtype.base, subshape)
        self.dtype = numpy.dtype([(n, d[n]) for n in d])

    def __contains__(self, column):
        return column in self.dtype.names

    def __getitem__(self, column):
        if column not in self.cache:
            self[column] = self.fetch(column) 
        return self.cache[column]

    def forget(self, column):
        del self.cache[column]
