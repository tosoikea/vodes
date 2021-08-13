

from functools import reduce
from typing import List
from vodes.symbolic.symbols import BoundedExpression, BoundedVariable
from vodes.symbolic.power import Power
from pymbolic.mapper import RecursiveMapper
from pymbolic.mapper.evaluator import UnknownVariableError
from pymbolic.primitives import Expression, Quotient, Power, Variable, Sum, Product, is_constant

class BoundedMapper(RecursiveMapper):
    def __init__(self, context: dict, symbol: BoundedVariable):
        assert(not context is None)
        assert(symbol)
        self.context = context
        self.symbol = symbol

    def __apply(self, op, l, r):
        bound = self.symbol.bound
        lexpr = l
        rexpr = r

        if isinstance(l, BoundedExpression):
            bound = bound.intersect(l.bound)
            lexpr = l.expr

        if isinstance(r, BoundedExpression):
            bound = bound.intersect(r.bound)
            rexpr = r.expr

        return op(lexpr, rexpr, bound)

    def _bop(self, l, r, b, op):
        r = op(l,r)

        if is_constant(r):
            return r
        else:
            return BoundedExpression(
                expression=r,
                boundary=b
            )

    def _badd(self, l, r, b):
        return self._bop(l,r,b,lambda a,b:Sum((a,b)))

    def _bmul(self, l, r, b):
        return self._bop(l,r,b,lambda a,b:Product((a,b)))

    def _bdiv(self, l, r, b):
        return self._bop(l,r,b,lambda a,b:Quotient(a,b))

    def _bpow(self, l, r, b):
        return self._bop(l,r,b,lambda a,b:Power(a,b))

    def map_constant(self, expr):
        return expr

    def map_variable(self, expr:Variable) -> List[BoundedExpression]:
        # we do not substitute the free symbol
        if not (self.symbol.name in self.context) and self.symbol.name == expr.name:
            return BoundedExpression(
                expr,
                boundary=self.symbol.bound
            )
        else:
            try:
                return self.context[expr.name]
            except KeyError:
                raise UnknownVariableError(expr.name)
    
    def map_sum(self, expr:Sum) -> BoundedExpression:
        return reduce(
            lambda r, x: self.__apply(self._badd, r, x), [
                self.rec(child) for child in expr.children
            ]
        )

    def map_product(self, expr:Product) -> Expression:
        return reduce(
            lambda r, x: self.__apply(self._bmul, r, x), [
                self.rec(child) for child in expr.children
            ]
        )

    def map_power(self, expr:Power) -> Expression:
        return self.__apply(self._bpow, self.rec(expr.base), self.rec(expr.exponent))

    def map_quotient(self, expr:Quotient):
        return self.__apply(self._bdiv, self.rec(expr.numerator), self.rec(expr.denominator))