from sys import intern

from vodes.symbolic.expressions.primitives import ExtendedExpression

# Expression Library

# Expression Mapper
from pymbolic.mapper.stringifier import StringifyMapper

class ConstantStringifyMapper(StringifyMapper):
    def map_pi(self, expr, enclosing_prec, *args, **kwargs):
        return f'π'

class Pi(ExtendedExpression):
    """The '\pi' (3.141592654\ldots) constant"""
    init_arg_names = ()

    def __init__(self):
        super().__init__()

    def __getinitargs__(self):
        return ()

    @property
    def expr(self):
        return self.expression

    def make_stringifier(self, originating_stringifier=None):
        return ConstantStringifyMapper()

    mapper_method = intern("map_pi")