import numpy

class Column(object):
    def __init__(self, parent, column):
        self.parent = parent
        self.column = column

    def __getitem__(self, index):
        if isinstance(index, slice):
            a, b, c = index.indices(self.parent.size)
            assert c == 1
            return self.parent._wrapped_fetch(self.column, a, b)
        else:
            a = int(index)
            b = a + 1
            return self.parent._wrapped_fetch(self.column, a, b)[0]

class Rows(object):
    def __init__(self, parent, rows):
        self.parent = parent
        self.rows = rows
    def __getitem__(self, index):
        return self.parent[index][self.rows]

class ColumnStore(object):
    """ A ColumnStore

        Columns are accessed via :py:code:`[columnname]`. 
        
        Subclass shall implement :py:meth:`fetch`, :py:attr:`dtype`, :py:attr:`size`.

        A ColumnStore acts as a context manager. During when the context is held,
        all fetch operations are buffered. The buffers are release after the context
        is released.

        Notes
        -----
        To retrive the contents of a columns use :py:code:`[columnname][:]`, or
        other slicing syntax with a step size of 1.

        >>> with mycolumnstore:
        >>>    print mycolumnstore['Column1'][:]
        >>>    print mycolumnstore['Column1'][:]
    """
    def __init__(self):
        self._memorybuffer_ = None

    @property
    def dtype(self):
        """dtype of each item"""
        raise NotImplementedError

    @property
    def size(self):
        """Total number of items """
        raise NotImplementedError

    def __enter__(self):
        self._memorybuffer_ = {}
        return self 

    def __exit__(self, a, b, c):
        self._memorybuffer_ = None

    def fetch(self, column, start, end):
        """Load data from start to end for a column """
        raise NotImplementedError

    def _wrapped_fetch(self, column, start, end):
        """Check _memorybuffer_ before fetch"""
        _memorybuffer_ = self._memorybuffer_
        if _memorybuffer_ is not None:
            if (column, start, end) in _memorybuffer_:
                data = _memorybuffer_[(column, start, end)]
            else:
                data = self.fetch(column, start, end)
                _memorybuffer_[(column, start, end)] = data
            return data
        else:
            return self.fetch(column, start, end)

    def __iter__(self):
        return iter(self.dtype.names)

    def __contains__(self, column):
        return column in self.dtype.names

    def __getitem__(self, column):
        if isinstance(column, slice):
            return Rows(self, column)
        if column in self:
            return Column(self, column)
        else:
            raise KeyError('column `%s` not found' % column)
