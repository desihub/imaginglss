from __future__ import print_function

__all__ = ['ColumnArray']

import numpy
import os
import os.path

class ColumnArray(object):
    """ a Column Array that can be stored as a 
        filesystem tree of binary files, plus some text format
        meta data.
    """
    def __init__(self, size=0):
        self.data = dict()
        self.size = size
        self.header = dict()
        
    def __len__(self):
        return self.size

    @property
    def dtype(self):
        items = [(k, self.data[k].dtype, self.data[k].shape) 
                    for k in self.data]
        items = [ (k, dtype if len(shape) == 1 else 
                  (dtype, shape[1:]))
            for k, dtype, shape in items]
        
        return numpy.dtype(items)

    def __setitem__(self, key, value):
        value = numpy.array(value, copy=False)
        if len(value) != self.size:
            raise ValueError("Size of new column does not match the columnArray")

        if key != key.strip() or ':' in key:
            raise KeyError("key must not contain ':' or white spaces")

        self.data[key] = value

    def __getitem__(self, index):
        """ for string index, returns a column.
            for iterable index, returns a ColumnArray viewing/copying
            only rows selected by the iterable.
        """
        if isinstance(index, basestring):
            return self.data[index]
        else:
            return self.select(index)

    def __contains__(self, key):
        return key in self.data

    def keys(self):
        return self.data.keys()

    def __iter__(self):
        return iter(self.data)

    def extend(self, newrows):
        """ grow the array with new rows """
        raise NotImplementedError

    def select(self, indices):
        """ returns a new ColumnArray with only items at indices.
        """
        if isinstance(indices, numpy.ndarray) and indices.dtype.char == '?':
            newsize = indices.sum()
        else:
            newsize = len(indices)
        subset = self.__class__(size=newsize)
        for key in self.data:
            subset[key] = self.data[key][indices]
        return subset 

    def as_ndarray(self):
        data = numpy.empty(self.size, dtype=self.dtype)
        for key in self:
            data[key][:] = self[key] 
        return data

    def __repr__(self):
        return "ColumnArray: columns =
    @classmethod
    def fromndarray(cls, array, copy=False, header=None):
        ca = cls(size=len(array))
        for key in array.dtype.names:
            ca[key] = array[key]
        return ca

    @classmethod
    def fromfile(cls, prefix):
        with file(os.path.join(prefix, '__dtype__.info')) as storage:
            lines = storage.readlines()

        assert lines[0].startswith('ROWS')
        size = int(lines[0].split(':')[1])
        ca = cls(size=size)
        for line in lines[1:]:
            key, dtype, shape = line.split(':')
            key = key.strip()
            dtype = numpy.dtype(dtype.strip())
            shape = parse_tuple(shape)
            ca[key] = numpy.fromfile(os.path.join(prefix, key), dtype=(dtype, shape))
        return ca

    def tofile(self, prefix):
        # first create the directory
        for key in self:
            filename = os.path.join(prefix, key)
            try:
                os.makedirs(os.path.dirname(filename))
            except OSError:
                pass
            self[key].tofile(filename)

        # now serialize the metadata in text
        with file(os.path.join(prefix, '__dtype__.info'), mode='w') as storage:
            dtype = self.dtype
            storage.write('ROWS : %d\n' % self.size)
            for key in self:
                storage.write('%s : %s : %s\n' % (key, dtype[key].base, str(dtype[key].shape)))

        with file(os.path.join(prefix, '__header__.info'), mode='w') as storage:
            header = self.header
            for key in header:
                storage.write('%s : %s\n' % (key, str(header[key])))

# SPAMS
import re
def parse_tuple(shape):
    """ parse (1, 2, 3) tuples.  """
    if not re.match('^[0-9, ()]*$', shape):
        raise ValueError('unsupported shape')
    return eval(shape)

def test():
    data1 = numpy.arange(10).reshape(5, 2)
    data2 = numpy.arange(100).reshape(5, 2, 10)
    data3 = numpy.arange(5)

    ca = ColumnArray(size=5)
    ca['col1d'] = data3
    ca['col2d'] = data1
    ca['col3d'] = data2
    try:
        ca['badcol'] = numpy.arange(3)
    except ValueError:
        pass
    assert 'col1d' in ca
    ca.tofile('testarray')

    ca2 = ColumnArray.fromfile('testarray')
    print(ca.dtype)
    print(ca2.dtype)
    print ca2

if __name__ == '__main__':
    test()
