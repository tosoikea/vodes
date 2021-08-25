

from functools import reduce
from typing import List
from vodes.symbolic.expressions.primitives import Subtraction

# Custom Expression Library
from vodes.symbolic.expressions.bounded import BoundedExpression, BoundedVariable
from vodes.symbolic.expressions.trigonometric import sin, cos

# Symbolic Expression Library
from pymbolic.mapper import RecursiveMapper
from pymbolic.mapper.evaluator import UnknownVariableError
from pymbolic.primitives import Expression, Quotient, Power, Variable, Sum, Product, is_constant

class BoundedMapper(RecursiveMapper):
    def __init__(self, context: dict, symbol: BoundedVariable):
        assert(not context is None)
        assert(symbol)
        self.context = context
        self.symbol = symbol

    def __apply(self, op, l, r=None):
        bound = self.symbol.bound
        lexpr = l
        rexpr = r

        if isinstance(l, BoundedExpression):
            bound = bound.intersect(l.bound)
            lexpr = l.expr

        if not (r is None) and isinstance(r, BoundedExpression):
            bound = bound.intersect(r.bound)
            rexpr = r.expr

        # unary
        if r is None:
            return op(lexpr, bound)
        else:
            return op(lexpr, rexpr, bound)

    def _bop(self, b, op, l, r=None):
        res = op(l) if r is None else op(l,r)

        if is_constant(res):
            return res
        else:
            return BoundedExpression(
                expression=res,
                boundary=b
            )

    def _badd(self, l, r, b):
        return self._bop(b,lambda a,b:Sum((a,b)),l,r=r)

    def _bsub(self, l, r, b):
        return self._bop(b,lambda a,b:Subtraction((a,b)),l,r=r)

    def _bmul(self, l, r, b):
        return self._bop(b,lambda a,b:Product((a,b)),l,r=r)

    def _bdiv(self, l, r, b):
        return self._bop(b,lambda a,b:Quotient(a,b),l,r=r)

    def _bpow(self, l, r, b):
        return self._bop(b,lambda a,b:Power(a,b),l,r=r)

    def _bsin(self, x, b):
        return self._bop(b,lambda a:sin(a),x)

    def _bcos(self, x, b):
        return self._bop(b,lambda a:cos(a),x)

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

    def map_sub(self, expr:Subtraction): 
        return reduce(
            lambda r, x: self.__apply(self._bsub, r, x), [
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

    def map_sin(self, expr):
        return self.__apply(self._bsin, self.rec(expr.expr))

    def map_cos(self, expr):
        return self.__apply(self._bcos, self.rec(expr.expr))