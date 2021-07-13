from sys import intern
from vodes.symbolic.symbols import Boundary, BoundedVariable
from pymbolic.mapper import RecursiveMapper
from pymbolic.mapper.stringifier import StringifyMapper
from pymbolic.primitives import Expression, Product, Sum, Variable


class Max(Expression):
    """Class to represent the maximum function as an expression.

    Args:
        expression: The expression to be wrapped within the maximum function.

    Attributes:
        expression: The wrapped expression.
    """
    init_arg_names = ("expression",)

    def __init__(self, expression):
        assert(not expression is None)

        self.expression = expression

    def __getinitargs__(self):
        return (self.expr,)

    @property
    def expr(self):
        return self.expression

    def make_stringifier(self, originating_stringifier=None):
        return MaxStringifyMapper()

    mapper_method = intern("map_maximum")

class MaxStringifyMapper(StringifyMapper):
    def map_maximum(self, expr, enclosing_prec, *args, **kwargs):
        return f'max({expr.expr})'