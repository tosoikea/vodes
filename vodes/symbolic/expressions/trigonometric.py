from sys import intern

from vodes.symbolic.expressions.primitives import ExtendedExpression

# Expression Library

# Expression Mapper
from pymbolic.mapper.stringifier import StringifyMapper

class TrigonometricStringifyMapper(StringifyMapper):
    def map_sin(self, expr, enclosing_prec, *args, **kwargs):
        return f'sin({expr.expr})'

    def map_cos(self, expr, enclosing_prec, *args, **kwargs):
        return f'cos({expr.expr})'

class sin(ExtendedExpression):
    """Class to represent the sinus function as an expression.

    Args:
        expression: The expression to be wrapped within the sinus function.

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
        return TrigonometricStringifyMapper()

    mapper_method = intern("map_sin")

class cos(ExtendedExpression):
    """Class to represent the sinus function as an expression.

    Args:
        expression: The expression to be wrapped within the sinus function.

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
        return TrigonometricStringifyMapper()

    mapper_method = intern("map_cos")
