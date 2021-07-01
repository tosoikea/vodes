from functools import reduce

from sys import intern
from pymbolic.primitives import Expression, Lookup, Variable
from pymbolic.mapper.evaluator import EvaluationMapper
from pymbolic.mapper.stringifier import PREC_NONE, StringifyMapper

##
# Auxiliary interval functions.
# src : https://link.springer.com/content/pdf/10.1007/3-540-36599-0_7.pdf
# src : https://www.math.kit.edu/ianm2/~kulisch/media/arjpkx.pdf
###
class Interval(Expression):
    init_arg_names = ("lower", "upper",)

    def __init__(self, lower, upper=None):
        # degenerate interval
        if upper is None:
            upper = lower

        self.lower = lower
        self.upper = upper

    def __getinitargs__(self):
        return self.lower, self.upper

    @property
    def low(self):
        return self.lower

    @property
    def up(self):
        return self.upper

    def make_stringifier(self, originating_stringifier=None):
        return IntervalStringifyMapper()

    mapper_method = intern("map_interval")

    # Inequality
    def __abs__(self):
        return max(abs(self.low), abs(self.up))

class IntervalEvaluator(EvaluationMapper):
    def __init__(self, context:dict, in_bounds=False):
        super().__init__(context=context)
        self.in_bounds=in_bounds

    def _iadd(self,l, r):
        return Interval(l.low + r.low, l.up + r.up)

    def _isub(self,l,r):
        return Interval(l.low - r.up, l.up - r.low)

    def _imul(self,l,r):
        lmin, lmax = l.lower, l.upper
        rmin, rmax = r.lower, r.upper
        
        # TODO : This can be simplified base on case distinction
        
        # Degenerate intervals
        if lmin == lmax:
            return Interval(lmin * rmin, lmin * rmax)
        elif rmin == rmax:
            return Interval(lmin * rmin, lmax * rmin)

        try:
            return Interval(
                min(lmin * rmin, lmin * rmax, lmax * rmin, lmax * rmax),
                max(lmin * rmin, lmin * rmax, lmax * rmin, lmax * rmax)
            )
        except TypeError:
            return l * r

    def _idiv(self,l,r):
        rmin, rmax = r.lower, r.upper

        # TODO : 0 \in Y assumed -> to be handled
        return self._imul(l,Interval(1/rmin, 1/rmax))

    #
    # TODO : Handle special cases (e.g. 0**-1)
    #
    def _ipow(self,l,r):
        lmin, lmax = l.lower, l.upper
        rmin, rmax = r.lower, r.upper

        try:
            return Interval(
                min(lmin ** rmin, lmin ** rmax, lmax ** rmin, lmax ** rmax),
                max(lmin ** rmin, lmin ** rmax, lmax ** rmin, lmax ** rmax)
            )
        except TypeError:
            return l ** r

    def _apply_operation(self, op, l, r):
        l = l if isinstance(l, Interval) else Interval(l)
        r = r if isinstance(r, Interval) else Interval(r)

        return op(l,r)

    def map_interval(self, expr):
        lower = self._eval_bound(expr.low)
        upper = self._eval_bound(expr.up)

        return Interval(lower, upper)

    def map_sum(self, expr):
        if not self.in_bounds:
            return self._map_iv_sum(expr)
        else:
            return super().map_sum(expr)

    def map_product(self, expr):
        if not self.in_bounds:
            return self._map_iv_product(expr)
        else:
            return super().map_product(expr)

    def map_quotient(self, expr):
        if not self.in_bounds:
            return self._map_iv_quotient(expr)
        else:
            return super().map_quotient(expr)
        
    # TODO : Only fabs supported
    def map_call(self, expr, *args):
        def make_f(name):
            return Lookup(Variable("math"), name)

        if expr.function == make_f('fabs'):
            children = [self.rec(par, *args) for i, par in enumerate(expr.parameters)]
            return abs(*children)
        else:
            raise NotImplementedError(f"The function {expr.function} is not supported")


    ##
    # Interval mapping. Can be overwritten, if desired.
    ##
    def _eval_bound(self, expr):
        return IntervalEvaluator(context=self.context, in_bounds=True)(expr)

    def _map_iv_product(self, expr):
        return reduce(lambda a,b: self._apply_operation(self._imul,a,b), [self.rec(child) for child in expr.children])

    def _map_iv_quotient(self, expr):
        return self._apply_operation(self._idiv, self.rec(expr.numerator), self.rec(expr.denominator))

    def _map_iv_sum(self, expr):
        return reduce(lambda a,b: self._apply_operation(self._iadd,a,b), [self.rec(child) for child in expr.children])
    # --


class IntervalStringifyMapper(StringifyMapper):
    def map_interval(self, expr, enclosing_prec, *args, **kwargs):
        lower = self.rec(expr.low, PREC_NONE, *args, **kwargs)
        upper = self.rec(expr.up, PREC_NONE, *args, **kwargs)

        return f'[{lower},{upper}]'
