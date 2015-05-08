import numpy

class Column(object):
    def __init__(self, parent, column):
        self.parent = parent
        self.column = column

    def __getitem__(self, index):
        if isinstance(index, slice):
            a, b, c = index.indices(self.parent.size)
            assert c == 1
            return self.parent.fetch(self.column, a, b)
        else:
            a = int(index)
            b = a + 1
            return self.parent.fetch(self.column, a, b)[0]

class ColumnStore(object):
    """ A cached column store 

        Subclass shall implement :py:method:`fetch` that loads data
        from a data source.

    """
    def __init__(self):
        pass

    @property
    def dtype(self):
        return None

    @property
    def size(self):
        return 0

    def __iter__(self):
        return iter(self.dtype.names)

    def __contains__(self, column):
        return column in self.dtype.names

    def __getitem__(self, column):
        if column in self:
            return Column(self, column)
        else:
            raise KeyError('column `%s` not found' % column)
