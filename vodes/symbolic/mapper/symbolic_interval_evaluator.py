import logging

from abc import abstractmethod, ABC
from typing import List

# Custom Expression library
from vodes.symbolic.expressions.bounded import BoundedVariable, BoundedExpression, Domain
from vodes.symbolic.expressions.interval import Interval
from vodes.symbolic.expressions.rational import Rational
from vodes.symbolic.expressions.nthroot import NthRoot
from vodes.symbolic.expressions.absolute import Abs

# Expression Library
from pymbolic.primitives import Quotient, Variable, Sum, Product, Expression, Power

# Expression Mapper
from pymbolic.mapper import RecursiveMapper

#TODO : Inward and outward rounding
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
        self._assumptions = {}

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

        self._logger.info(f'(CONST) : {expr} -> {list(map(str,res))}')
        return res

    def map_rational(self, expr:Rational) -> List[BoundedExpression]:
        # Rational = Constant num and den
        # However, pymbolic has bad support for rationals -> Convert
        res = [
            BoundedExpression(
                Interval(
                    Quotient(
                        expr.numerator,
                        expr.denominator
                    )
                ),
                boundary=self._symbol.bound
            )
        ]
        
        self._logger.info(f'(RATIONAL) : {expr} -> {list(map(str,res))}')
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

        self._logger.info(f'(VAR) : {expr} -> {list(map(str,res))}')
        return res

    def map_sum(self, expr:Sum) -> List[BoundedExpression]:
        from functools import reduce
        res = reduce(
            lambda r, x: self.__apply(self._iadd, r, x), [
                self.rec(child) for child in expr.children
            ]
        )

        self._logger.info(f'(SUM) : {expr} -> {list(map(str,res))}')
        return res

    def map_product(self, expr:Product) -> List[BoundedExpression]:
        from functools import reduce
        children = [
                self.rec(child) for child in expr.children
            ]

        res = reduce(
            lambda r, x: self.__apply(self._imul, r, x), children
        )

        self._logger.info(f'(PRODUCT) : {expr} -> {list(map(str,res))}')
        return res

    def map_quotient(self, expr:Quotient) -> List[BoundedExpression]:
        res = self.__apply(
            self._idiv,
            self.rec(expr.numerator),
            self.rec(expr.denominator)
        )

        self._logger.info(f'(QUOTIENT) : {expr} -> {list(map(str,res))}')
        return res

    def map_power(self, expr:Power) -> List[BoundedExpression]:
        res = self.__apply(
            self._ipow,
            self.rec(expr.base),
            self.rec(expr.exponent)
        )

        self._logger.info(f'(POWER) : {expr} -> {list(map(str,res))}')
        return res
    
    def map_nthroot(self, expr:NthRoot) -> List[BoundedExpression]:
        bexprs = self.rec(expr.expr)
        res = []

        for bexpr in bexprs:
            res.extend(
                self. _inthroot(
                    bexpr.expr,
                    n=expr.n,
                    d=bexpr.bound
                )
            )

        self._logger.info(f'(NTHROOT) : {expr} -> {list(map(str,res))}')
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

        self._logger.info(f'(INTERVAL) {expr} -> {list(map(str,res))}')
        return res

    def map_absolute(self, expr:Abs) -> List[BoundedExpression]:
        res = []
        bexprs = self.rec(expr.expr)

        for bexpr in bexprs:
            self._logger.warn(bexpr)
            res.extend(
                self._iabs(i=bexpr.expr, d=bexpr.bound)
            )

        self._logger.info(f'(ABS) : {expr} -> {list(map(str,res))}')
        return res

    ####
    # INTERVAL INTERFACE
    ####
    @abstractmethod
    def _iadd(self, l:Interval, r:Interval, d: Domain) -> List[BoundedExpression]:
        """Interval addition as defined within https://epubs.siam.org/doi/abs/10.1137/1.9780898717716.ch2.
        
        Args:
            l (Interval): The left parameter of the addition.
            r (Interval): The right parameter of the addition.
            b (Boundary): The boundary ('Gültigkeitsbereich') of the addition.

        Returns:
            _iadd: A list of BoundedExpressions containing the result of the addition (interval), possibly limited to subsets of the initial boundary, that may in sum not be equal to the initial boundary.
        """
        pass

    @abstractmethod
    def _isub(self, l:Interval, r:Interval, d: Domain) -> List[BoundedExpression]:
        """Interval substitution as defined within https://epubs.siam.org/doi/abs/10.1137/1.9780898717716.ch2.
        
        Args:
            l (Interval): The left parameter of the substitution.
            r (Interval): The right parameter of the substitution.
            b (Boundary): The boundary ('Gültigkeitsbereich') of the substitution.

        Returns:
            _isub: A list of BoundedExpressions containing the result of the substitution (interval), possibly limited to subsets of the initial boundary, that may in sum not be equal to the initial boundary.
        """
        pass

    @abstractmethod
    def _idiv(self, l:Interval, r:Interval, d: Domain) -> List[BoundedExpression]:
        """Interval division as defined within https://epubs.siam.org/doi/abs/10.1137/1.9780898717716.ch2.
        
        Args:
            l (Interval): The dividend of the division.
            r (Interval): The divisor of the division.
            b (Boundary): The boundary ('Gültigkeitsbereich') of the division.

        Returns:
            _idiv: A list of BoundedExpressions containing the result of the division (interval), possibly limited to subsets of the initial boundary, that may in sum not be equal to the initial boundary.
        """
        pass

    @abstractmethod
    def _imul(self, l:Interval, r:Interval, d: Domain) -> List[BoundedExpression]:
        """Interval multiplication as defined within https://epubs.siam.org/doi/abs/10.1137/1.9780898717716.ch2. Because symbolic intervals are supported, the inheriting classes define custom evaluations for symbolic cases.
        
        Args:
            l (Interval): The left parameter of the multiplication.
            r (Interval): The right parameter of the multiplication.
            b (Boundary): The boundary ('Gültigkeitsbereich') of the multiplication.

        Returns:
            _imul: A list of BoundedExpressions containing the result of the multiplication (interval), possibly limited to subsets of the initial boundary.
        """
        pass

    @abstractmethod
    def _ipow(self, l:Interval, r:Interval, d: Domain) -> List[BoundedExpression]:
        pass

    @abstractmethod
    def _inthroot(self, i:Interval, n:int, d:Domain) -> List[BoundedExpression]:
        pass
    
    @abstractmethod
    def _iabs(self, i:Interval, d:Domain) -> List[BoundedExpression]:
        pass
    

    @abstractmethod
    def _icontains(self, expr, val, incl, bf, d: Domain) -> list:
        """Determine the inclusion of any values of a boundary based on the inclusion of the upper and lower boundaries. For the interval to include a value, both boundaries have to include it."""  
        pass

    ####
    # SYMBOLIC EXPRESSION INTERFACE
    ####
    @abstractmethod
    def _minimum(self, exprs:List[Expression],boundary:Domain) -> List[BoundedExpression]:
        pass

    @abstractmethod
    def _maximum(self, exprs:List[Expression],boundary:Domain) -> List[BoundedExpression]:
        pass

    ####
    # UTILITY FUNCTIONS
    ####
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

    def contains(self, iv:Interval, xs, d: Domain) -> list:
        """Determine the inclusion of any values of a boundary based on the inclusion of the upper and lower boundaries. For the interval to include a value, both boundaries have to include it."""  
        from vodes.symbolic.utils import le,ge
        from vodes.symbolic.expressions.infinity import NegativeInfinity, Infinity

        if xs is None:
            return [
                (
                    d,
                    []
                )
            ]

        if not isinstance(xs, Domain):
            bs = Domain(xs,xs,False,False)
        else:
            bs = xs

        lower_sections = []
        upper_sections = []

        (lv,lo) = (bs.start,bs.left_open) if not isinstance(bs.start, NegativeInfinity) else (bs.end,bs.right_open)
        (rv,ro) = (bs.end,bs.right_open) if not isinstance(bs.end, Infinity) else (bs.start,bs.left_open)

        if not isinstance(bs.start, NegativeInfinity):
            self._logger.debug(f'Analyzing upper interval boundary {iv.up} for inclusion of {rv}')
            upper_sections = [
                (b,[xs] if len(vals) > 0 else []) for (b,vals) in self._icontains(
                    expr=iv.up,
                    val=rv,
                    incl=lambda y, val : ge(y,val),
                    bf=lambda included: (ro and included) or (not ro and not included),
                    d=d
                )
            ]
        else:
            upper_sections = [(d,[xs])]

        
        if not isinstance(bs.end, Infinity):
            self._logger.debug(f'Analyzing lower interval boundary {iv.low} for inclusion of {lv}')
            lower_sections = [
                (b,[xs] if len(vals) > 0 else []) for (b,vals) in self._icontains(
                    expr=iv.low,
                    val=lv,
                    incl=lambda y, val : le(y,val),
                    # if bf = true -> interval open
                    bf=lambda included: (lo and included) or (not lo and not included),
                    d=d
                )
            ]
        else:
            lower_sections = [(d,[xs])]

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