from sys import intern

# Custom Expression Library
from vodes.symbolic.expressions.bounded import Domain, BoundedVariable

# Expression Library
from pymbolic.primitives import Expression, Variable

# Expression Mapper
from pymbolic.mapper import RecursiveMapper
from pymbolic.mapper.stringifier import StringifyMapper


class NthRoot(Expression):
    """Class to represent the nth root function as an expression.

    Args:
        expression: The expression to be wrapped within the absolute function.

    Attributes:
        expression: The wrapped expression.
    """
    init_arg_names = ("n","expression",)

    def __init__(self, expression, n:int):
        assert(not expression is None)
        assert n > 1

        self._expression = expression
        self._n = n

    def __getinitargs__(self):
        return (self.expr,self.n)

    @property
    def expr(self):
        return self._expression

    @property
    def n(self):
        return self._n

    def make_stringifier(self, originating_stringifier=None):
        return NthRootStringifyMapper()

    mapper_method = intern("map_nthroot")

class NthRootStringifyMapper(StringifyMapper):
    def map_nthroot(self, expr, enclosing_prec, *args, **kwargs):
        return f'{expr.n}^√({expr.expr})'