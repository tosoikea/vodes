import logging
from typing import List, Tuple
from vodes.symbolic.expressions.infinity import NegativeInfinity
from vodes.symbolic.translations.nop_translation import NOPTranslation

# Assumption library
from vodes.symbolic.translations.to_scalar import ToScalar
from vodes.symbolic.properties.is_scalar import IsScalar
from vodes.symbolic.assumption import Assumption

# Custom Expression library
from vodes.symbolic.expressions.interval import Interval
from vodes.symbolic.expressions.bounded import BoundedVariable, BoundedExpression, Domain
from vodes.symbolic.expressions.rational import Rational
from vodes.symbolic.expressions.primitives import Subtraction
from vodes.symbolic.expressions.absolute import Abs
from vodes.symbolic.expressions.nthroot import NthRoot
from vodes.symbolic.expressions.trigonometric import sin, cos

# Custom Mappers
from vodes.symbolic.mapper.simplification_mapper import simplify
from vodes.symbolic.mapper.interop import ExactPymbolicToSympyMapper, ExactSympyToPymbolicMapper

# Expression Library
from pymbolic.primitives import Expression, Quotient, Variable, Sum, Product, Power

# Expression Mapper
from pymbolic.mapper import RecursiveMapper

def evaluate(expression, context=None, float:bool=False):
    if context is None:
        context = {}
    res = ScalarEvaluator(context)(expression)

    # Two iterations of solver, if symbolic values are used for evaluation.
    # This allows to push the floating calculations further up.
    if float:
        return evaluate(res, context=context)
    else:
        return res

class ScalarEvaluator(RecursiveMapper):
    """Class for determining the exact boundaries of intervals on the basis of function analysis."""
    @classmethod
    def is_multivariant(cls) -> bool:
        return True

    def __init__(self, context: dict):
        assert(not context is None)

        self._logger = logging.getLogger(__name__)
        self._context = context

    ####
    # EXPRESSION MAPPING
    ####
    def map_constant(self, expr) -> Expression:
        res = Interval(expr)
        self._logger.debug(f'(CONST) : {expr} -> {res}')
        return res
  
    def map_variable(self, expr:Variable) -> Expression:
        from vodes.symbolic.mapper.extended_evaluation_mapper import evaluate
        from pymbolic.mapper.evaluator import UnknownVariableError

        res = None

        try:
            subs = self._context[expr.name]

            if isinstance(subs,Interval):
                res = subs
            else:
                res = Interval(subs)
        except KeyError:
            raise UnknownVariableError(expr.name)

        self._logger.debug(f'(VAR) : {expr} -> {res}')
        return res

    def map_sum(self, expr:Sum) -> Expression:
        from functools import reduce
        res = ExactSympyToPymbolicMapper()(
            ExactPymbolicToSympyMapper()(
                reduce(lambda r,x: self._iadd(r,x),[self.rec(child) for child in expr.children])
            )
        )

        self._logger.debug(f'(SUM) : {expr} -> {res}')
        return res

    def map_sub(self, expr:Subtraction) -> Expression:
        from functools import reduce
        res = reduce(lambda r,x: self._isub(r,x),[self.rec(child) for child in expr.children])

        self._logger.debug(f'(SUB) : {expr} -> {res}')
        return res

    def map_product(self, expr:Product) -> Expression:
        from functools import reduce
        res = reduce(lambda r,x: self._imul(r,x),[self.rec(child) for child in expr.children])

        self._logger.debug(f'(PROD) : {expr} -> {res}')
        return res

    def map_quotient(self, expr:Quotient) -> Expression:
        from vodes.symbolic.mapper.interop import ExactPymbolicToSympyMapper, ExactSympyToPymbolicMapper
        
        res = self._idiv(
            self.rec(expr.numerator),
            self.rec(expr.denominator)
        )
        self._logger.debug(f'(SUB) : {expr} -> {res}')
        return res

    def map_power(self, expr:Power) -> Expression:
        res = self._ipow(
            self.rec(expr.base),
            self.rec(expr.exponent)
        )

        self._logger.debug(f'(POW) : {expr} -> {res}')
        return res

    def map_interval(self, expr:Interval) -> Expression:
        l = self.rec(expr.low)
        r = self.rec(expr.up)

        exprs = [
            l.low,
            l.up,
            r.low,
            r.up
        ]

        res = Interval(
            self._minimum(exprs),
            self._maximum(exprs)
        )

        self._logger.debug(f'(INT) : {expr} -> {res}')
        return res

    def map_absolute(self, expr:Abs) -> Expression:
        res = self._iabs(self.rec(expr.expr))

        self._logger.debug(f'(ABS) : {expr} -> {res}')
        return res

    ## FUNCTIONS
    def map_nthroot(self, expr:NthRoot) -> Expression:
        res = self._inthroot(self.rec(expr.expr),n=expr.n)
        
        self._logger.debug(f'(NTHROOT) : {expr} -> {res}')
        return res

    def map_sin(self, expr:sin) -> Expression:
        res = self._isin(self.rec(expr.expr))
        
        self._logger.debug(f'(SIN) : {expr} -> {res}')
        return res

    def map_cos(self, expr:cos) -> Expression:
        res = self._icos(self.rec(expr.expr))
        
        self._logger.debug(f'(COS) : {expr} -> {res}')
        return res

    ####
    # INTERVAL INTERFACE
    ####
    def _iadd(self, l:Interval, r:Interval) -> Expression:
        """Interval addition as defined within https://epubs.siam.org/doi/abs/10.1137/1.9780898717716.ch2.
        
        Args:
            l (Interval): The left parameter of the addition.
            r (Interval): The right parameter of the addition.
            b (Boundary): The boundary ('Gültigkeitsbereich') of the addition.

        Returns:
            _iadd: A list of BoundedExpressions containing the result of the addition (interval), possibly limited to subsets of the initial boundary, that may in sum not be equal to the initial boundary.
        """
        return Interval(l.low + r.low, l.up + r.up)

    def _isub(self, l:Interval, r:Interval) -> Expression:
        """Interval substitution as defined within https://epubs.siam.org/doi/abs/10.1137/1.9780898717716.ch2.
        
        Args:
            l (Interval): The left parameter of the substitution.
            r (Interval): The right parameter of the substitution.
            b (Boundary): The boundary ('Gültigkeitsbereich') of the substitution.

        Returns:
            _isub: A list of BoundedExpressions containing the result of the substitution (interval), possibly limited to subsets of the initial boundary, that may in sum not be equal to the initial boundary.
        """
        return Interval(l.low - r.up, l.up - r.low)

    def _idiv(self, l:Interval, r:Interval) -> Expression:
        """Interval division as defined within https://epubs.siam.org/doi/abs/10.1137/1.9780898717716.ch2.
        
        Args:
            l (Interval): The dividend of the division.
            r (Interval): The divisor of the division.
            b (Boundary): The boundary ('Gültigkeitsbereich') of the division.

        Returns:
            _idiv: A list of BoundedExpressions containing the result of the division (interval), possibly limited to subsets of the initial boundary, that may in sum not be equal to the initial boundary.
        """
        from vodes.symbolic.utils import le,ge
        rmin, rmax = r.low, r.up

        # 1. 0 \in r
        if le(rmin, 0) and ge(rmax, 0):
            # Infinity = Not Defined (for us)
            self._logger.warning(f'Dropping {r} for division, as it contains a 0.')
            return None
        # 2. 0 \not\in r
        else:
            return self._imul(l, Interval(Quotient(1,rmax),Quotient(1,rmin)))

    def _imul(self, l:Interval, r:Interval) -> Expression:
        """Interval multiplication as defined within https://epubs.siam.org/doi/abs/10.1137/1.9780898717716.ch2. Because symbolic intervals are supported, the inheriting classes define custom evaluations for symbolic cases.
        
        Args:
            l (Interval): The left parameter of the multiplication.
            r (Interval): The right parameter of the multiplication.
            b (Boundary): The boundary ('Gültigkeitsbereich') of the multiplication.

        Returns:
            _imul: A list of BoundedExpressions containing the result of the multiplication (interval), possibly limited to subsets of the initial boundary.
        """
        exprs = [
                l.low * r.low, 
                l.low * r.up, 
                l.up * r.low, 
                l.up * r.up
            ]

        lower = self._minimum(exprs)
        upper = self._maximum(exprs)

        return Interval(
            lower,
            upper
        )

    def _ipow(self, l:Interval, r:Interval) -> Expression:
        from vodes.symbolic.mapper.interop import ExactPymbolicToSympyMapper, ExactSympyToPymbolicMapper

        if not (r.low) == (r.up):
            raise ValueError("Not supporting non degenerate exponents")

        return ExactSympyToPymbolicMapper()(
            ExactPymbolicToSympyMapper()(l) ** ExactPymbolicToSympyMapper()(r)
        )
    
    def _iabs(self, i:Interval) -> Expression:
        return Interval(self._maximum(
            exprs=[
                i.up,
                (-1) * i.low
            ]
        ))

    def _inthroot(self, i:Interval, n:int) -> Expression:
        # linear
        return Interval(
            NthRoot(i.low,n),
            NthRoot(i.up,n)
        )

    def _isin(self, i:Interval) -> Expression:
        pass

    def _icos(self, i:Interval) -> Expression:
        pass

    ####
    # SYMBOLIC EXPRESSION INTERFACE
    ####
    def _minimum(self, exprs:List[Expression]) -> Expression:
        from vodes.symbolic.utils import minimum
        res = minimum(exprs)

        self._logger.debug(f'Minimum {res} determined from {list(map(str,exprs))}')
        return res

    def _maximum(self, exprs:List[Expression]) -> Expression:
        from vodes.symbolic.utils import maximum

        res = maximum(exprs)

        self._logger.debug(f'Maximum {res} determined from {list(map(str,exprs))}')
        return res

    