import logging

from abc import abstractmethod, ABC
from typing import Dict, List
from vodes.symbolic.expressions.primitives import Subtraction

# Custom Expression library
from vodes.symbolic.expressions.bounded import BoundedVariable, BoundedExpression, Domain
from vodes.symbolic.expressions.interval import Interval
from vodes.symbolic.expressions.rational import Rational
from vodes.symbolic.expressions.nthroot import NthRoot
from vodes.symbolic.expressions.trigonometric import cos, sin
from vodes.symbolic.expressions.absolute import Abs

# Expression Library
from pymbolic.primitives import Quotient, Variable, Sum, Product, Expression, Power

# Expression Mapper
from pymbolic.mapper import RecursiveMapper

class IntervalEvaluator(ABC, RecursiveMapper):
    """Class to evaluate expressions containing symbolic intervals. These expressions have to be limited to one free variable after substitution.

    Args:
        context (dict): The values to substitute for symbols within the expressions. Has to include all free variables, except the variable specified using symbol.
        symbol (BoundedVariable): The variable used for the symbolic intervals. Its value has to be bounded.

    Attributes:
        _context (dict): The stored substitution dictionary.
        _symbol (BoundedVariable): The only allowed free variables.
    """
    def __init__(self, context: dict, symbol: BoundedVariable):
        assert (symbol)

        self._logger = logging.getLogger(__name__)

        self._context = context if context is not None else {}
        self._symbol = symbol
        self._assumptions = {
            "_iadd" : ([],[]),      #l,r
            "_isub" : ([],[]),      #l,r
            "_idiv" : ([],[]),      #l,r
            "_imul" : ([],[]),      #l,r
            "_ipow" : ([],[]),      #l,r
            "_iabs" : ([]),         #l
            "_inthroot" : ([]),     #l
            "_isin" : ([]),         #l
            "_icos" : ([]),         #l
            "_icontains" : ([]),    #l
            "_minimum" : ([]),      #l
            "_maximum" : ([])       #l
        }

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
        from vodes.symbolic.mapper.extended_evaluation_mapper import evaluate
        from pymbolic.mapper.evaluator import UnknownVariableError

        res = None

        # vexpr do not substitute free symbols
        if not expr.name in self._context and expr.name == self._symbol.name:
            vexpr = Interval(expr)
        else:
            try:
                subs = self._context[expr.name]

                if isinstance(subs,Interval):
                    vexpr = subs
                else:
                    vexpr = Interval(subs)
            except KeyError:
                raise UnknownVariableError(expr.name)

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

    def map_sub(self, expr:Subtraction) -> List[BoundedExpression]:
        from functools import reduce
        res = reduce(
            lambda r, x: self.__apply(self._isub, r, x), [
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

    def map_interval(self, expr:Interval) -> List[BoundedExpression]:
        from vodes.symbolic.utils import merge

        lower = self.rec(expr.low)
        upper = self.rec(expr.up)

        res = [
            BoundedExpression(
                expression=Interval(
                    lower=left.low,
                    upper=right.up
                ),
                boundary=boundary
            ) for ((left,right),boundary) in merge(lower,upper)
        ]

        self._logger.info(f'(INTERVAL) {expr} -> {list(map(str,res))}')
        return res

    def map_absolute(self, expr:Abs) -> List[BoundedExpression]:
        res = self.__unary_apply(
            self._iabs,
            self.rec(expr.expr)
        )

        self._logger.info(f'(ABS) : {expr} -> {list(map(str,res))}')
        return res

    ## FUNCTIONS
    def map_nthroot(self, expr:NthRoot) -> List[BoundedExpression]:
        res = self.__unary_apply(
            self._inthroot,
            self.rec(expr.expr),
            opt=(expr.n,)
        )

        self._logger.info(f'(NTHROOT) : {expr} -> {list(map(str,res))}')
        return res

    def map_sin(self, expr:sin) -> List[BoundedExpression]:
        raise NotImplementedError()

    def map_cos(self, expr:cos) -> List[BoundedExpression]:
        raise NotImplementedError()

    ####
    # INTERVAL INTERFACE
    ####
    @abstractmethod
    def _iadd(self, l:Interval, r:Interval, d:Domain) -> List[BoundedExpression]:
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
    def _isub(self, l:Interval, r:Interval, d:Domain) -> List[BoundedExpression]:
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
    def _idiv(self, l:Interval, r:Interval, d:Domain) -> List[BoundedExpression]:
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
    def _imul(self, l:Interval, r:Interval, d:Domain) -> List[BoundedExpression]:
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
    def _ipow(self, l:Interval, r:Interval, d:Domain) -> List[BoundedExpression]:
        pass
    
    @abstractmethod
    def _iabs(self, i:Interval, d:Domain) -> List[BoundedExpression]:
        pass

    @abstractmethod
    def _inthroot(self, i:Interval, n:int, d:Domain) -> List[BoundedExpression]:
        pass

    @abstractmethod
    def _isin(self, i:Interval, d:Domain) -> List[BoundedExpression]:
        pass

    @abstractmethod
    def _icos(self, i:Interval, d:Domain) -> List[BoundedExpression]:
        pass
    
    @abstractmethod
    def _icontains(self, expr:Interval, val, d:Domain, incl:set=set(("up","low"))) -> list:
        """Determine the inclusion of any values of a boundary based on the inclusion of the upper and lower boundaries. For the interval to include a value, both boundaries have to include it."""  
        pass

    ####
    # SYMBOLIC EXPRESSION INTERFACE
    ####
    @abstractmethod
    def _minimum(self, exprs:List[Expression],d:Domain) -> List[BoundedExpression]:
        pass

    @abstractmethod
    def _maximum(self, exprs:List[Expression],d:Domain) -> List[BoundedExpression]:
        pass

    ####
    # UTILITY FUNCTIONS
    ####
    def _apply_assumptions(self, assumptions:list, exprs:list, d):
        translated = [BoundedExpression(expression=expr,boundary=d) for expr in exprs]

        for assumption in assumptions:
            translated = assumption.validate(
                translated
            )

            # TODO
            # At the moment, we expect the expression at the index i of exprs to equal the expression at the same index after the translation.
            # If we were to support division within the translations, we require lists of lists to not loose the association between the initial and the translated expression(s).
            # We therefore introduce the assumption : exprs[i] => res[i].expr && boundary == res[i].bound
            #  
            if len(translated) != len(exprs):
                raise ValueError("It is not yet supported, to divide the domain of expressions within translations.")
        
        return [val.expr for val in translated]

    def __convert_expr(self, l, d):
        expr = l

        if isinstance(l, BoundedExpression):
            d = d.intersect(l.bound)
            expr = l.expr
        elif isinstance(l, Variable) or self.is_constant(l):
            expr = Interval(l)

        # Implicit interval conversion (degenerate)
        if not isinstance(expr, Interval):
            expr = Interval(expr)

        return (expr,d)

    def __unary_apply(self, op, ls:list, opt:tuple = None):
        res = []        
        
        (lassums) = self._assumptions.get(op.__name__)

        for l in ls:
            (lexpr,domain) = self.__convert_expr(l,self._symbol.bound)

            if domain is None:
                continue

            lexpr = self._apply_assumptions(lassums,[lexpr],domain)[0]

            res.extend(
                op(lexpr, domain) if opt is None else op(lexpr, *opt, domain)
            )

        return res

    def __apply(self, op, ls: list, rs: list):
        res = []

        (lassums,rassums) = self._assumptions.get(op.__name__)

        for l in ls:
            for r in rs:
                (lexpr,domain) = self.__convert_expr(l,self._symbol.bound)
                (rexpr,domain) = self.__convert_expr(r,domain)

                if domain is None:
                    continue

                res.extend(
                    op(
                        self._apply_assumptions(lassums,[lexpr],domain)[0], 
                        self._apply_assumptions(rassums,[rexpr],domain)[0],
                        domain
                    )
                )

        return res

    def contains(self, iv:Interval, xs, d:Domain) -> list:
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

        (lv,lo) = (bs.start,bs.left_open) if not isinstance(bs.start, NegativeInfinity) else (bs.end,bs.right_open)
        (rv,ro) = (bs.end,bs.right_open) if not isinstance(bs.end, Infinity) else (bs.start,bs.left_open)

        if not isinstance(bs.start, NegativeInfinity):
            self._logger.debug(f'Analyzing upper interval boundary {iv.up} for inclusion of {rv}')
            upper_sections = [
                (b,[xs] if len(vals) > 0 else []) for (b,vals) in self._icontains(
                    expr=iv,
                    val=rv,
                    incl=set(("up",)),
                    d=d
                )
            ]
        else:
            upper_sections = [(d,[xs])]

        if not isinstance(bs.end, Infinity):
            self._logger.debug(f'Analyzing lower interval boundary {iv.low} for inclusion of {lv}')
            lower_sections = [
                (b,[xs] if len(vals) > 0 else []) for (b,vals) in self._icontains(
                    expr=iv,
                    val=lv,
                    incl=set(("low",)),
                    d=d
                )
            ]
        else:
            lower_sections = [(d,[xs])]

        res = []
        
        # Correct inclusions
        # if bf = true -> interval open
        bf_lower=lambda included: (lo and included) or (not lo and not included)
        bf_upper=lambda included: (ro and included) or (not ro and not included) 

        for i in range(len(lower_sections)):
            (lower_b, lower_vals) = lower_sections[i]

            llo = bf_lower(xs in lower_vals) if i > 0 else d.left_open
            lro = bf_lower(xs in lower_vals) if i < len(lower_sections) - 1 else d.right_open
            
            for j in range(len(upper_sections)):
                (upper_b, upper_vals) = upper_sections[j]

                rlo = bf_upper(xs in upper_vals) if j > 0 else d.left_open
                rro = bf_upper(xs in upper_vals) if j < len(upper_sections) - 1 else d.right_open

                lower_b = Domain(
                    start=lower_b.start,
                    end=lower_b.end,
                    left_open=llo,
                    right_open=lro)

                upper_b = Domain(
                    start=upper_b.start,
                    end=upper_b.end,
                    left_open=rlo,
                    right_open=rro)

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