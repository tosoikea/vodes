from sys import intern
from pymbolic.mapper.stringifier import StringifyMapper
from pymbolic.primitives import Expression


class Rational(Expression):
    """Class to represent a rational number.

    Args:
        numerator: The positive integer numerator
        denominator: The positive integer denominator

    Attributes:
        num: Numerator of the rational number.
        den: Denominator of the rational number.
    """
    init_arg_names = ("numerator","denominator",)

    def __init__(self, numerator, denominator):
        assert(isinstance(numerator,int) and numerator >= 0)
        assert(isinstance(denominator,int) and denominator > 0)

        self.numerator = numerator
        self.denominator = denominator

    def __getinitargs__(self):
        return (self.expr,)

    @property
    def num(self):
        return self.numerator

    @property
    def den(self):
        return self.denominator

    def make_stringifier(self, originating_stringifier=None):
        return RationalStringifyMapper()

    mapper_method = intern("map_rational")

class RationalStringifyMapper(StringifyMapper):
    def map_absolute(self, expr, enclosing_prec, *args, **kwargs):
        return f'{expr.num} / {expr.den}'