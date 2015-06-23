"""
npyquery, Query mini-language for numpy arrays

Examples
--------
>>> d = dtype([
        ('BlackholeMass', 'f4'), 
        ('Position', ('f4', 3))])

>>> data = numpy.zeros(4, dtype=d)
>>> data['BlackholeMass'][:] = numpy.arange(4)
>>> data['Position'][:] = numpy.arange(12).reshape(3, 1)
>>> query  = (Column('BlackholeMass') > 2.0)
>>> query &= (Column('BlackholeMass') > 4.0)
>>> query &= (Column('Position')[:, 2] > 1.0)
>>> print query.visit(data)

"""
import numpy

class Node(object):
    def __le__(self, other):
        return Expr("<=", numpy.less_equal, [self, other])
    def __lt__(self, other):
        return Expr("<", numpy.less, [self, other])
    def __eq__(self, other):
        return Expr("==", numpy.equal, [self, other])
    def __gt__(self, other):
        return Expr(">", numpy.greater, [self, other])
    def __ge__(self, other):
        return Expr(">=", numpy.greater_equal, [self, other])
    def __and__(self, other):
        return Expr("&", numpy.bitwise_and, [self, other])
    def __or__(self, other):
        return Expr("|", numpy.bitwise_or, [self, other])
    def __xor__(self, other):
        return Expr("^", numpy.bitwise_xor, [self, other])
    def __invert__(self):
        return Expr("~", numpy.bitwise_not, [self])
    def __pow__(self, other):
        return Expr("**", numpy.power, [self, other])
    def __rpow__(self, other):
        return Expr("**", numpy.power, [other, self])
    def __mul__(self, other):
        return Expr("*", numpy.multiply, [self, other])
    def __div__(self, other):
        return Expr("/", numpy.divide, [self, other])
    def __add__(self, other):
        return Expr("+", numpy.add, [self, other])
    def __neg__(self):
        return Expr("-", numpy.negative, [self])
    def __sub__(self, other):
        return Expr("-", numpy.subtract, [self, other])
    def __mod__(self, other):
        return Expr("%", numpy.remainder, [self, other])
    def __getitem__(self, index):
        return GetItem(self, index)
    def sin(self):
        return Expr("sin", numpy.sin, [self])
    def cos(self):
        return Expr("cos", numpy.cos, [self])
    def tan(self):
        return Expr("tan", numpy.tan, [self])
    def log(self):
        return Expr("log", numpy.log, [self])
    def log10(self):
        return Expr("log10", numpy.log10, [self])
    def __call__(self, array):
        return array[self.visit(array)]
    def visit(self, array):
        raise NotImplemented 

def repr_slice(s):
    if not isinstance(s, slice):
        return repr(s)
    if s.start is None and s.stop is None and s.step is None:
        return ':'
    start = "" if s.start is None else repr(s.start)
    stop = "" if s.stop is None else repr(s.stop)
    if s.step is None:
        return "%s:%s" % (start, stop)
    else:
        step = repr(s.step)
    return "%s:%s:%s" % (start, stop, step)

class GetItem(Node):
    def __init__(self, obj, index):
        self.obj = obj
        self.index = index
    def __repr__(self):
        if isinstance(self.index, tuple):
            rslice = ','.join([repr_slice(o) for o in self.index])
        else:
            rslice = repr_slice(self.index)
        return "%s[%s]" % (repr(self.obj), rslice)

    def visit(self, array):
        return self.obj.visit(array)[self.index]

class Literal(Node):
    def __init__(self, value):
        self.value = value
    def __repr__(self):
        return repr(self.value)
    def visit(self, array):
        return self.value

class Column(Node):
    def __init__(self, column):
        Node.__init__(self)
        self.column = column
    def __repr__(self):
        return "[%s]" % self.column
    def visit(self, array):
        return array[self.column]

class Expr(Node):
    def __init__(self, operator, function, operands):
        self.operator = operator
        self.function = function
        operands = [
            a if isinstance(a, Node)
            else Literal(a)
            for a in operands
        ]

        self.operands = operands
        self.flatten()
        
    def is_associative(self):
        if not isinstance(self.function, numpy.ufunc):
            return False
        if self.function.identity is not None:
            return True
        return False

    def __iter__(self):
        return iter(self.operands)

    def flatten(self):
        if not self.is_associative(): return
        o = []
        for a in self.operands:
            if not isinstance(a, Expr):
                o.append(a)
                continue
            if a.function != self.function:
                o.append(a)
                continue
            else:
                o.extend(a.operands)
        self.operands = o

    def __repr__(self):
        if len(self.operands) >= 2:
            return "(%s)" % (' ' + self.operator + ' ').join([repr(a) for a in self.operands])
        elif len(self.operands) == 1:
            [a] = self.operands
            return "(%s %s)" % (self.operator, repr(a))
        else:
            raise ValueError

    def visit(self, array):
        ops = [a.visit(array) for a in self.operands]
        if self.is_associative():
            return self.function.reduce(ops)
        else:
            return self.function(*ops)

def test():    
    d = numpy.dtype([
        ('BlackholeMass', 'f4'), 
        ('PhaseOfMoon', 'f4'), 
        ('Position', ('f4', 3)),
])

    data = numpy.zeros(5, dtype=d)
    data['BlackholeMass'][:] = numpy.arange(len(data))
    data['PhaseOfMoon'][:] = numpy.linspace(0, 1, len(data), endpoint=True)
    data['Position'][:] = numpy.arange(data['Position'].size).reshape(len(data), -1)
    query  = (Column('BlackholeMass') > 0.0)
    query &= (Column('BlackholeMass') < 5.0)
    query &= (Column('Position')[:, 2] > 0.0) | (Column('Position')[:, 1] < 0.0)
    query &= (numpy.sin(Column('PhaseOfMoon') * (2 * numpy.pi)) < 0.1)
    print query
    print query(data)
    for sub in query:
        print sub, sub.visit(data)
if __name__ == '__main__':
    test()