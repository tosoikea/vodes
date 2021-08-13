import logging

from typing import List
from abc import abstractmethod, ABC

from sys import intern

from vodes.symbolic.absolute import Abs
from vodes.symbolic.maximum import Max
from vodes.symbolic.power import Power
from vodes.symbolic.rational import Rational
from vodes.symbolic.symbols import Boundary, BoundedExpression, BoundedValue, BoundedVariable, Infinity, NegativeInfinity

from pymbolic.mapper import RecursiveMapper
from pymbolic.primitives import Expression, Quotient, Variable, Sum, Product
from pymbolic.mapper.stringifier import PREC_NONE, StringifyMapper


##
# src : https://link.springer.com/content/pdf/10.1007/3-540-36599-0_7.pdf
# src : https://www.math.kit.edu/ianm2/~kulisch/media/arjpkx.pdf
###
class Interval(Expression):
    """Class to represent an interval, as defined within interval analysis, as expression.

    Args:
        lower: The lower boundary of the interval, may be symbolic.
        upper: The upper boundary of the interval. If none is given, a degenerate interval [lower,lower] is constructed.

    Attributes:
        __lower: The lower boundary of the interval. Use property low!
        __upper: The upper boundary of the interval. Use property up!
    """
    init_arg_names = ("lower", "upper",)

    def __init__(self, lower, upper=None):
        assert(not lower is None)

        # degenerate interval
        if upper is None:
            upper = lower

        self.__lower = lower
        self.__upper = upper

    def __getinitargs__(self):
        return self.low, self.up

    @property
    def low(self):
        """Get or set the lower boundary of the interval"""
        return self.__lower

    @property
    def up(self):
        """Get or set the upper boundary of the interval"""
        return self.__upper

    def make_stringifier(self, originating_stringifier=None):
        return IntervalStringifyMapper()

    mapper_method = intern("map_interval")

class IntervalStringifyMapper(StringifyMapper):
    def map_interval(self, expr, enclosing_prec, *args, **kwargs):
        lower = self.rec(expr.low, PREC_NONE, *args, **kwargs)
        upper = self.rec(expr.up, PREC_NONE, *args, **kwargs)

        return f'[{lower},{upper}]'


class SymbolicIntervalEvaluator(ABC, RecursiveMapper):
    """Class to evaluate expressions containing symbolic intervals. These expressions have to be limited to one free variable after substitution.

    Args:
        context (dict): The values to substitute for symbols within the expressions. Has to include all free variables, except the variable specified using symbol.
        symbol (BoundedVariable): The variable used for the symbolic intervals. Its value has to be bounded.

    Attributes:
        _context (dict): The stored substitution dictionary.
        _symbol (BoundedVariable): The stored symbol, only allowed free variable.
    """
    def __init__(self, context: dict, symbol: BoundedVariable):
        assert(not context is None)
        assert(symbol)

        self._logger = logging.getLogger(__name__)
        self._context = context
        self._symbol = symbol

    ####
    # EXPRESSION MAPPING
    ####
    def map_constant(self, expr) -> List[BoundedExpression]:
        res = [
            BoundedExpression(
                Interval(expr),
                boundary=self._symbol.bound
            )
        ]

        self._logger.debug(f'{expr} -> {list(map(str,res))}')
        return res

    def map_rational(self, expr:Rational) -> List[BoundedExpression]:
        # Rational = Constant num and den
        # However, pymbolic has bad support for rationals -> Convert
        res = [
            BoundedExpression(
                Interval(
                    expr.num * Power(expr.den,-1)
                ),
                boundary=self._symbol.bound
            )
        ]
        
        self._logger.debug(f'{expr} -> {list(map(str,res))}')
        return res
  
    def map_variable(self, expr:Variable) -> List[BoundedExpression]:
        from pymbolic import evaluate

        res = None
        # vexpr do not substitute the free symbol
        if not str(self._symbol.name) in self._context and self._symbol.name == expr.name:
            vexpr = Interval(expr)
        else:
            vexpr = Interval(evaluate(expression=expr,context=self._context))

        res = [
            BoundedExpression(
                expression=vexpr,
                boundary=self._symbol.bound
            )
        ]

        self._logger.debug(f'{expr} -> {list(map(str,res))}')
        return res

    def map_sum(self, expr:Sum) -> List[BoundedExpression]:
        from functools import reduce
        res = reduce(
            lambda r, x: self.__apply(self._iadd, r, x), [
                self.rec(child) for child in expr.children
            ]
        )

        self._logger.debug(f'{expr} -> {list(map(str,res))}')
        return res

    def map_product(self, expr:Product) -> List[BoundedExpression]:
        from functools import reduce
        res = reduce(
            lambda r, x: self.__apply(self._imul, r, x), [
                self.rec(child) for child in expr.children
            ]
        )

        self._logger.debug(f'{expr} -> {list(map(str,res))}')
        return res

    def map_quotient(self, expr:Quotient) -> List[BoundedExpression]:
        res = self.__apply(
            self._idiv,
            self.rec(expr.numerator),
            self.rec(expr.denominator)
        )

        self._logger.debug(f'{expr} -> {list(map(str,res))}')
        return res

    def map_power(self, expr:Power) -> List[BoundedExpression]:
        res = self.__apply(
            self._ipow,
            self.rec(expr.base),
            self.rec(expr.exponent)
        )

        self._logger.debug(f'{expr} -> {list(map(str,res))}')
        return res

    def map_interval(self, expr:Interval) -> List[BoundedExpression]:
        from vodes.symbolic.mapper.bounded_mapper import BoundedMapper
        l = BoundedMapper(context=self._context, symbol=self._symbol)(expr.low)
        r = BoundedMapper(context=self._context, symbol=self._symbol)(expr.up)

        lexpr = l
        rexpr = r

        if isinstance(l, BoundedExpression):
            lexpr = l.expr

        if isinstance(r, BoundedExpression):
            rexpr = r.expr

        res = [
            BoundedExpression(
                expression=Interval(lower=lexpr, upper=rexpr),
                boundary=self._symbol.bound
            )
        ]
        self._logger.debug(f'{expr} -> {list(map(str,res))}')
        return res

    def map_maximum(self, expr:Max) -> List[BoundedExpression]:
        bexprs = self.rec(expr.expr)

        res = [
            BoundedExpression(
                expression=item.expr.up,
                boundary=item.boundary
            ) 
            for item in bexprs
        ]
        self._logger.debug(f'{expr} -> {list(map(str,res))}')
        return res

    def map_absolute(self, expr:Abs) -> List[BoundedExpression]:
        bexprs = self.rec(expr.expr)
        res = [
            item for sublist in list(map(lambda bexpr : self._sabs(bexpr.expr, bexpr.bound), bexprs))
                 for item in sublist
        ]

        self._logger.debug(f'{expr} -> {list(map(str,res))}')
        return res

    ####
    # INTERVAL INTERFACE
    ####
    def _iadd(self, l, r, b: Boundary):
        """Interval addition as defined within https://epubs.siam.org/doi/abs/10.1137/1.9780898717716.ch2.
        
        Args:
            l (Interval): The left parameter of the addition.
            r (Interval): The right parameter of the addition.
            b (Boundary): The boundary ('Gültigkeitsbereich') of the addition.

        Returns:
            _iadd: A single element list, containing the result of the addition (interval) within an BoundedExpression.
        """
        return [
            BoundedExpression(
                expression=Interval(l.low + r.low, l.up + r.up),
                boundary=b
            )
        ]

    def _isub(self, l, r, b: Boundary):
        """Interval substitution as defined within https://epubs.siam.org/doi/abs/10.1137/1.9780898717716.ch2.
        
        Args:
            l (Interval): The left parameter of the substitution.
            r (Interval): The right parameter of the substitution.
            b (Boundary): The boundary ('Gültigkeitsbereich') of the substitution.

        Returns:
            _isub: A single element list, containing the result of the substitution (interval) within an BoundedExpression.
        """
        return [
            BoundedExpression(
                expression=Interval(l.low - r.up, l.up - r.low),
                boundary=b
            )
        ]

    def _idiv(self, l, r, b: Boundary):
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

        if self.is_constant(rmin) and self.is_constant(rmax):
            # 1. 0 \in r
            if le(rmin, 0) and ge(rmax, 0):
                # Infinity = Not Defined (for us)
                self._logger.debug(f'Dropping {r} for division, as it contains a 0.')
                return []
            # 2. 0 \not\in r
            else:
                return self._imul(l, Interval(Quotient(1,rmax),Quotient(1,rmin)), b)
        else:
            return self._sdiv(l, r, b)

    def _imul(self, l, r, b: Boundary):
        """Interval multiplication as defined within https://epubs.siam.org/doi/abs/10.1137/1.9780898717716.ch2. Because symbolic intervals are supported, the inheriting classes define custom evaluations for symbolic cases.
        
        Args:
            l (Interval): The left parameter of the multiplication.
            r (Interval): The right parameter of the multiplication.
            b (Boundary): The boundary ('Gültigkeitsbereich') of the multiplication.

        Returns:
            _imul: A list of BoundedExpressions containing the result of the multiplication (interval), possibly limited to subsets of the initial boundary.
        """
        res = None
        lmin, lmax = l.low, l.up
        rmin, rmax = r.low, r.up

        # TODO : Simplify assertion
        if self.is_constant(lmin) and self.is_constant(lmax) and self.is_constant(rmin) and self.is_constant(rmax):
            exprs = [
                lmin * rmin, 
                lmin * rmax, 
                lmax * rmin, 
                lmax * rmax
            ]

            (ilow,iup) = self.__constant_min_max(exprs)

            res = Interval(
                ilow,
                iup
            )
        else:
            # min, max could not be handled => result of symbolic boundaries
            return self._smul(l, r, b)

        return [
            BoundedExpression(
                expression=res,
                boundary=b
            )
        ]
    
    def _icontains(self, iv:Interval, xs, b: Boundary) -> list:
        """Determine the inclusion of any values of a boundary based on the inclusion of the upper and lower boundaries. For the interval to include a value, both boundaries have to include it."""  

        lower_sections = []
        upper_sections = []


        if not isinstance(xs, Boundary):
            bs = Boundary(BoundedValue(xs,False),BoundedValue(xs,False))
        else:
            bs = xs


        lower_val = bs.lower if not isinstance(bs.lower.value, NegativeInfinity) else bs.upper
        upper_val = bs.upper if not isinstance(bs.upper.value, Infinity) else bs.lower

        if not isinstance(bs.lower.value, NegativeInfinity):
            self._logger.debug(f'Analyzing upper interval boundary {iv.up} for inclusion of {upper_val}')
            upper_sections = [
                (b,[xs] if len(vals) > 0 else []) for (b,vals) in self._sinclusion(
                    iv.up,
                    upper_val.value,
                    incl=lambda y, val : y >= val,
                    b=b,
                    bf=lambda included: (upper_val.open and included) or (not upper_val.open and not included)
                )
            ]
        else:
            upper_sections = [(b,[xs])]

        
        if not isinstance(bs.upper.value, Infinity):
            self._logger.debug(f'Analyzing lower interval boundary {iv.low} for inclusion of {lower_val}')
            lower_sections = [
                (b,[xs] if len(vals) > 0 else []) for (b,vals) in self._sinclusion(
                    iv.low,
                    lower_val.value,
                    incl=lambda y, val : y <= val,
                    b=b,
                    # if bf = true -> interval open
                    bf=lambda included: (lower_val.open and included) or (not lower_val.open and not included)
                )
            ]
        else:
            lower_sections = [(b,[xs])]

        res = []
        for (lower_b, lower_vals) in lower_sections:
            for (upper_b, upper_vals) in upper_sections:
                combined = lower_b.intersect(upper_b)

                if combined is None:
                    continue

                # One includes, other excludes => In summary excluded
                included = [v for v in lower_vals if v in upper_vals]

                res.append(
                    (
                        combined,
                        included
                    )
                )
            
        self._logger.debug(f'Combined boundary section {lower_sections}(l) and {upper_sections}(u) to {res}')
        return res


    ## TODO : Evaluate concept
    def _ipow(self, l, r, b: Boundary):
        return self._spow(l, r, b)

    ####
    # SYMBOLIC INTERVAL INTERFACE
    ####
    @abstractmethod
    def _smul(self, l:Interval, r:Interval, b: Boundary):
        """Interval multiplication in the case of symbolic interval boundaries. In this case, the normal min/max interval multiplication can not be evaluated.

        Args:
            l (Interval): The left parameter of the multiplication.
            r (Interval): The right parameter of the multiplication.
            b (Boundary): The boundary ('Gültigkeitsbereich') of the multiplication.

        Returns:
            _smul: A list of BoundedExpressions containing the result of the symbolic multiplication (interval), possibly limited to subsets of the initial boundary.
        """
        pass

    @abstractmethod
    def _sdiv(self, l:Interval, r:Interval, b: Boundary):
        """Interval division in the case of symbolic interval boundaries. In this case, the normal interval divison can not be determined, as knowledge about the inclusion of 0 is required.

        Args:
            l (Interval): The dividend of the division.
            r (Interval): The divisor of the division.
            b (Boundary): The boundary ('Gültigkeitsbereich') of the multiplication.

        Returns:
            _sdiv: A list of BoundedExpressions containing the result of the symbolic division (interval), possibly limited to subsets of the initial boundary.
        """
        pass

    @abstractmethod
    def _spow(self, l:Interval, r:Interval, b: Boundary):
        pass

    @abstractmethod
    def _sabs(self, i:Interval, b: Boundary):
        pass

    @abstractmethod
    def _sinclusion(self, expr, val, incl, bf, b: Boundary) -> list:
        pass

    ####
    # UTILITY FUNCTIONS
    ####
    def is_constant(self, expr):
        """Determines, if the expression is constant"""
        from pymbolic.primitives import is_constant
        from pymbolic.mapper.dependency import DependencyMapper
        return is_constant(expr) or not bool(DependencyMapper()(expr))

    def __apply(self, op, ls: list, rs: list):
        res = []

        for l in ls:
            for r in rs:
                bound = self._symbol.bound
                lexpr = l
                rexpr = r

                if isinstance(l, BoundedExpression):
                    bound = bound.intersect(l.bound)
                    lexpr = l.expr
                # May result from BoundedMapper
                elif isinstance(l, Variable) or self.is_constant(l):
                    lexpr = Interval(l)

                if isinstance(r, BoundedExpression):
                    bound = bound.intersect(r.bound)
                    rexpr = r.expr
                # May result from BoundedMapper
                elif isinstance(r, Variable) or self.is_constant(l):
                    rexpr = Interval(l)

                if not (isinstance(lexpr, Interval) and isinstance(rexpr, Interval)):
                    raise TypeError(
                        f"Cannot process {lexpr}, {rexpr} as (atleast) one entity is not an interval.")

                res.extend(op(lexpr, rexpr, bound))

        return res


    def __constant_min_max(self, exprs:list) -> tuple:
        """Determine the minimum and maximum for constant pymbolic expressions"""
        from pymbolic import evaluate
        from pymbolic.primitives import is_constant
        
        r_min = exprs[0]
        r_max = exprs[0]

        for expr in exprs:
            if is_constant(expr):
                lt = expr < r_min
            else:
                lt = evaluate(expr.lt(r_min))

            if lt:
                r_min = expr
                continue

            if is_constant(expr):
                gt = expr > r_max
            else:
                gt = evaluate(expr.lt(r_max))

            if gt:
                r_max = expr

        return (r_min,r_max)


class ApproxIntervalEvaluator(SymbolicIntervalEvaluator):
    """Super class for those interval evaluators, who evaluate symbolic cases approximatively.
    """
    def __init__(self, context: dict, symbol: Variable):
        super().__init__(context=context, symbol=symbol)


class ExactIntervalEvaluator(SymbolicIntervalEvaluator):
    """Super class for those interval evaluators, who evaluate symbolic cases exactly (costly).
    """
    def __init__(self, context: dict, symbol: Variable):
        super().__init__(context=context, symbol=symbol)