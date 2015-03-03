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
        self.updatedtype(column, value)

    def updatedtype(self, column, value):
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
            value = self.fetch(column) 
            self[column] = value.astype(self.dtype[column])
        return self.cache[column]

    def __delitem__(self, column):
        del self.cache[column]

import os.path
class DiskColumnStore(ColumnStore):
    def __init__(self, root, dtype=None):
        self.root = root
        ColumnStore.__init__(self, dtype)

    def getfilename(self, column):
        return os.path.join(self.root, column)
    
    def __getitem__(self, column):
        if column not in self.cache:
            filename = self.getfilename(column)
            try:
                print 'reading'
                data = numpy.fromfile(filename,
                    dtype=self.dtype[column])
                print data.dtype
                print 'reading done'
            except IOError:
                print 'fetching'
                data = ColumnStore.__getitem__(self, column)
                try:
                    os.makedirs(os.path.dirname(filename))
                except OSError:
                    pass
                print 'writing'
                print data.dtype
                data.tofile(filename)
                print 'writing done'
            self[column] = data
        return self.cache[column]

    def forget(self, column):
        try:
            os.unlink(self.getfilename(column))
        except IOError:
            pass
        del self[column]


    
