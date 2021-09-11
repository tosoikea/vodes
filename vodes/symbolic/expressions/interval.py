from sys import intern
from pymbolic.mapper.stringifier import PREC_NONE, StringifyMapper
from vodes.symbolic.expressions.primitives import ExtendedExpression

##
# src : https://link.springer.com/content/pdf/10.1007/3-540-36599-0_7.pdf
# src : https://www.math.kit.edu/ianm2/~kulisch/media/arjpkx.pdf
###
class Interval(ExtendedExpression):
    """Class to represent an interval, as defined within interval analysis, as expression.

    Args:
        lower: The lower boundary of the interval, may be symbolic.
        upper: The upper boundary of the interval. If none is given, a degenerate interval [lower,lower] is constructed.

    Attributes:
        __lower: The lower boundary of the interval. Use property low!
        __upper: The upper boundary of the interval. Use property up!
    """
    init_arg_names = ("_lower", "_upper",)

    def __init__(self, lower, upper=None):
        assert(not lower is None)

        # degenerate interval
        if upper is None:
            upper = lower

        self._lower = lower
        self._upper = upper

    def __getinitargs__(self):
        return self.low, self.up

    @property
    def low(self):
        """Get or set the lower boundary of the interval"""
        return self._lower

    @property
    def up(self):
        """Get or set the upper boundary of the interval"""
        return self._upper

    def make_stringifier(self, originating_stringifier=None):
        return IntervalStringifyMapper()

    mapper_method = intern("map_interval")

class IntervalStringifyMapper(StringifyMapper):
    def map_interval(self, expr, enclosing_prec, *args, **kwargs):
        lower = self.rec(expr.low, PREC_NONE, *args, **kwargs)
        upper = self.rec(expr.up, PREC_NONE, *args, **kwargs)

        return f'[{lower},{upper}]'