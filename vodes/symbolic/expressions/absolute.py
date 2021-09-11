from sys import intern

# Custom Expression Library
from vodes.symbolic.expressions.bounded import Domain, BoundedVariable
from vodes.symbolic.expressions.primitives import ExtendedExpression

# Expression Library
from pymbolic.primitives import Variable

# Expression Mapper
from pymbolic.mapper import RecursiveMapper
from pymbolic.mapper.stringifier import StringifyMapper


class Abs(ExtendedExpression):
    """Class to represent the absolute function as an expression.

    Args:
        expression: The expression to be wrapped within the absolute function.

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
        return AbsoluteStringifyMapper()

    mapper_method = intern("map_absolute")

class AbsoluteStringifyMapper(StringifyMapper):
    def map_absolute(self, expr, enclosing_prec, *args, **kwargs):
        return f'|{expr.expr}|'