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
    """ A ColumnStore

        Columns are accessed via :py:code:`[columnname]`. 
        
        Subclass shall implement :py:meth:`fetch`, :py:attr:`dtype`, :py:attr:`size`.

        Notes
        -----
        To retrive the contents of a columns use :py:code:`[columnname][:]`, or
        other slicing syntax with a step size of 1.

    """
    def __init__(self):
        pass

    @property
    def dtype(self):
        """dtype of each item"""
        raise NotImplementedError

    @property
    def size(self):
        """Total number of items """
        raise NotImplementedError

    def fetch(self, column, start, end):
        """Load data from start to end for a column """
        raise NotImplementedError
    
    def __iter__(self):
        return iter(self.dtype.names)

    def __contains__(self, column):
        return column in self.dtype.names

    def __getitem__(self, column):
        if column in self:
            return Column(self, column)
        else:
            raise KeyError('column `%s` not found' % column)
