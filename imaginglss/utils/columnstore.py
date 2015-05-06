import numpy

class ColumnStore(object):
    """ A cached column store 

        Subclass shall implement :py:method:`fetch` that loads data
        from a data source.

    """
    def __init__(self, dtype=None):
        self.cache = {}
        if dtype is not None:
            self.dtype = dtype

    def __setitem__(self, column, value):
        value = numpy.array(value)
        if column in self.dtype.fields:
            self.cache[column] = value.astype(self.dtype[column].base)
        else:
            self.updatedtype(column, value)

    def updatedtype(self, column, value):
        d = dict(self.dtype.fields)
        subshape = value.shape[1:]
        if len(subshape) == 0:
            d[column] = (value.dtype.base, ())
        else:
            d[column] = (value.dtype.base, subshape)
        self.dtype = numpy.dtype([(n, d[n]) for n in d])

    def __iter__(self):
        return iter(self.dtype.names)

    def __contains__(self, column):
        return column in self.dtype.names

    def __getitem__(self, column):
        if column not in self.cache:
            value = self.fetch(column) 
            self[column] = value
        return self.cache[column]

    def __delitem__(self, column):
        del self.cache[column]

import os.path
class DiskColumnStore(ColumnStore):
    def __init__(self, root, dtype):
        self.root = root
        ColumnStore.__init__(self, dtype)

    def getfilename(self, column):
        return os.path.join(self.root, column)
    
    def __getitem__(self, column):
        if column not in self.cache:
            filename = self.getfilename(column)
            try:
                data = numpy.fromfile(filename,
                    dtype=self.dtype[column])
            except IOError:
                data = ColumnStore.__getitem__(self, column)
                try:
                    os.makedirs(os.path.dirname(filename))
                except OSError:
                    pass
                data.tofile(filename)
            self[column] = data
        return self.cache[column]

    def forget(self, column):
        try:
            os.unlink(self.getfilename(column))
        except OSError:
            pass
        try:
            del self[column]
        except KeyError:
            pass

