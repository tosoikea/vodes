import logging
from typing import List

# Custom Expression
from vodes.symbolic.expressions.interval import Interval
from vodes.symbolic.expressions.bounded import BoundedExpression, Domain
from vodes.symbolic.expressions.infinity import Infinity, NegativeInfinity

# Custom Mapper
from vodes.symbolic.mapper.interop import ExactPymbolicToSympyMapper

# Symbolic Expression Library
from pymbolic.primitives import Expression

# Mapper
from pymbolic.interop.sympy import SympyToPymbolicMapper

class Analysis:
    def __init__(self, pf:Expression, d:Domain=None) -> None:
        from vodes.symbolic.mapper.comparison_evaluator import ComparisonEvaluator as Evaluator

        assert not (pf is None)
        self.__pf = pf
        self.__f = ExactPymbolicToSympyMapper()(self.__pf)

        # TODO : Complex support
        self.__domain = Domain(NegativeInfinity(),Infinity(),left_open=True,right_open=True) if d is None else d

        if len(self.__f.free_symbols) > 1:
            raise ValueError("Not supporting functions with more than one free variable. Consider encoding the parameter as a vector.")

        if len(self.__f.free_symbols) == 1:
            self.__symbol = list(self.__f.free_symbols)[0]

            # We use an interval evaluator, to determine proper symbolic interval constructions.
            # TODO : Make this configurable
            self._evaluator = Evaluator(context={},symbol=self.symbol)
        else:
            self.__symbol = None
            self._evaluator = None

        self._logger = logging.getLogger(__name__)
        self.__limit = 10**3

    @property
    def func(self):
        """Get the problem of this analysis"""
        return self.__f

    @property
    def expr(self):
        """Get the problem of this analysis"""
        return self.__pf

    @property
    def symbol(self):
        """Get the free variable of the problem"""
        return self.__symbol

    def equals(self, y, d:Domain=None) -> List[Domain]:
        ask = self.__f - y

        if d is None:
            d =  self.__domain

        return self._solve(ask,d=d)

    def compare(self, problem, comparison, d:Domain=None):
        """Compare this problem with the expression y using the supplied comparison."""
        if isinstance(problem,Analysis):
            cmp_prob = problem
        else:
            cmp_prob = Analysis(problem)

        if d is None:
            d =  self.__domain

        ask = comparison(self.func,cmp_prob.func)
        return self._solve(ask,d=d)

    def evaluate(self, x):
        """Evaluate a function at a given point"""

        if isinstance(x,Expression):
            x = ExactPymbolicToSympyMapper()(x)

        # ~ is_costant() of sympy objects
        if self.__symbol is None:
            return self.func

        return self.__f.subs(self.__symbol, x)

    def limit(self, x):
        """Evaluate the limit of a function"""
        from sympy.series.limits import limit

        # ~ is_costant() of sympy objects
        if self.__symbol is None:
            return self.func

        return limit(self.__f, self.__symbol, x)

    def diff(self,n:int=1):
        from sympy import diff

        # ~ is_costant() of sympy objects
        if self.__symbol is None:
            return self.func

        res = self.func

        for i in range(n):
            res = diff(res,self.__symbol)

        return res

    def __abs_maximum(self,f):
        from sympy import diff

        # Boundary values
        vals = [
            self.__domain.start,
            self.__domain.end
        ]

        f_diff = diff(f,self.symbol)

        for r in self._solve(f_diff,self.__domain):
            if r.start != r.end:
                self._logger.warning(f'Dropping {r} as domains are not yet fully supported')

            vals.append(r.start)

        maximum = None
        for v in vals:
            res = abs(f.subs(self.symbol,v))
            maximum = res if maximum is None else max(maximum,res)
        
        return maximum

    def taylor(self, a=0, n:int=1) -> List[BoundedExpression]:
        """Approximates the given problem using a taylor expansion around the supplied a and the idea of taylor models.
        Given the problem f(x) this function returns a symbolic interval I(x) with the property \forall_x f(x) \in I(x).
        
        Args:
            a: Center of the taylor expansion. Defaults to 0.
            n: Degree of the taylor polynomial. Defaults to 1.

        Raises:
            TODO : Not (n+1) times differentiable on domain, f^(n) not continuous on domain

        Returns:
            taylor: Symbolic interval I(x) with taylor polynomials of degree n as boundary.
        """
        from vodes.symbolic.utils import merge
        from sympy import series, lambdify, Pow, Float
        from scipy.optimize import minimize_scalar
        from mpmath import mpf
        from vodes.symbolic.mapper.extended_evaluation_mapper import evaluate
        from math import factorial
        from fractions import Fraction

        assert n > 0

        if isinstance(a,Expression):
            a = ExactPymbolicToSympyMapper()(a)

        # To reduce abs(x-a)^(n+1) to (x-a)^(n+1)
        if n % 2 != 1 and a != 0:
            raise ValueError(f"The combination of exponent {n} with expansion point {a} is not yet supported.")

        # Constant Expression -> Exactly representable by order 1 or even 0 taylor polynomial
        if self.symbol is None:
            return [
                BoundedExpression(
                    expression=Interval(self.func),
                    boundary=self.__domain
                )
            ]

        # I(x) = [T_n(x),T_n(x)] + [-C,+C]
        # Use the taylor remainder estimation theorem :
        # C    = M / (n+1)! * |x-a|^(n+1)
        # M   >= | R_n(x) |, \forall_x 

        # (1) Construct the taylor polynomial 
        t_n = series(self.func,self.symbol,x0=a,n=n+1).removeO()

        # (2) Bound the remainder R_n(x)
        r_n = self.diff(n=n+1)

        # Constant expression
        if len(r_n.free_symbols) == 0:
            m = abs(r_n)
        else:

            # Analytical approach
            try:
                m = self.__abs_maximum(r_n)
            except ValueError:
                self._logger.warning(f'Falling back to optimizers, as {r_n} could not be solved analytically.')

                # we want to maximize
                minimize_func = -abs(r_n)

                r_n_func = lambdify(list(r_n.free_symbols),minimize_func)
                if isinstance(self.__domain.start,NegativeInfinity) and isinstance(self.__domain.end,Infinity):
                    res = minimize_scalar(r_n_func)
                else:
                    # TODO : We may limit our boundary here using the evaluation
                    b_left = None if isinstance(self.__domain.start, NegativeInfinity) else evaluate(self.__domain.start)
                    b_right = None if isinstance(self.__domain.end, Infinity) else evaluate(self.__domain.end)
                    res = minimize_scalar(r_n_func, bounds=(b_left, b_right), method='bounded')

                if not res.success:
                    raise ValueError("Could not determine optimum.")

                m = r_n_func(res.x)
                if m == mpf("inf") or m == mpf("-inf"):
                    raise ValueError("Could not reliably bound the error. Try to use the domain to restrict the functional values.")

        # TODO : Under & Overflow -> Determine the error accordingly
        m = float(m)
        m_approx = Fraction(m).limit_denominator(self.__limit)

        # This is allowed, as we only increase the size of our constant M while still upholding the inequality M >= |R_n(x)|
        if abs(m) < (1/self.__limit) and m != 0:
            m_approx = Fraction(1,self.__limit)

        #TODO Verallgemeinern für (n+1) uneven
        c = m_approx * Pow(factorial(n+1),-1) * (self.symbol - a)**(n+1)

        # (3) Construct valid symbolic intervals
        exprs = [
            SympyToPymbolicMapper()(t_n - c), 
            SympyToPymbolicMapper()(t_n + c)
        ]

        lower = self._evaluator._minimum(exprs,self.__domain)
        upper = self._evaluator._maximum(exprs,self.__domain)

        return [
            BoundedExpression(
                expression=Interval(
                    lower=left,
                    upper=right
                ),
                boundary=boundary
            ) for (left,right,boundary) in merge(lower,upper)
        ]

    def is_polynomial(self,n:int):
        from sympy import Poly
        from sympy.polys.polyerrors import PolynomialError

        # ~ is_costant() of sympy objects
        if self.__symbol is None:
            return True

        try:
            res = Poly(self.func, self.__symbol)
            return res.degree() <= n
        except PolynomialError:
            return False

    ####
    # UTILITY FUNCTIONS
    ####
    def _convert_solution(self,d:Domain):
        from vodes.symbolic.expressions.constants import Pi
        from sympy.core.numbers import Pi as SymPi, NegativeInfinity as SymNegInf, Infinity as SymInf
        from sympy.core.basic import Basic
        from sympy import CRootOf
        
        if isinstance(d.start,CRootOf):
            raise ValueError(f'Encountered unsolved root {d.start}')
        
        if isinstance(d.end,CRootOf):
            raise ValueError(f'Encountered unsolved root {d.start}')

        if isinstance(d.start,Basic):
            # Constant mapping, TODO : Custom SympyToPymbolicMapper
            if isinstance(d.start,SymPi):
                d.start = Pi()
            # Constant mapping, TODO : Custom SympyToPymbolicMapper
            elif isinstance(d.start,SymNegInf):
                d.start = NegativeInfinity()
            # Constant mapping, TODO : Custom SympyToPymbolicMapper
            elif isinstance(d.start,SymInf):
                d.start = Infinity()
            else:
                d.start = SympyToPymbolicMapper()(d.start)
        
        if isinstance(d.end,Basic):
            # Constant mapping, TODO : Custom SympyToPymbolicMapper
            if isinstance(d.end,SymPi):
                d.end = Pi()
            # Constant mapping, TODO : Custom SympyToPymbolicMapper
            elif isinstance(d.end,SymNegInf):
                d.end = NegativeInfinity()
            # Constant mapping, TODO : Custom SympyToPymbolicMapper
            elif isinstance(d.end,SymInf):
                d.end = Infinity()
            else:
                d.end = SympyToPymbolicMapper()(d.end)
        
    def _purge_solutions(self, ds:List[Domain], d:Domain=None) -> list:
        from sympy import im
        res = []

        for _d in ds:
            #TODO : Make real/complex support a parameter choice
            #TODO : Define this within the domains itself -> They then specify the support for complex
            #As a workaround we check if boundaries are imaginary
            if im(_d.end) != 0 or im(_d.start) != 0:
                self._logger.debug(f'Dropping {_d}.')
                continue

            self._convert_solution(_d)

            if (d is None):
                res.append(_d)
            else:
                combined = _d.intersect(d)
                if combined is None:
                    continue
                else:
                    res.append(combined)

        return res

    def _map_solution(self, xs, d:Domain=None) -> list:
        from sympy.sets.sets import Union, Intersection, Complement, Interval as SymInterval
        from sympy.sets.fancysets import Reals

        res = []

        if xs.is_empty:
            self._logger.debug(f'No solution.')
        elif xs.is_FiniteSet:
            res = []
            self._logger.debug(xs)
            res = self._purge_solutions(ds=[
                Domain(start=x) for x in xs
            ], d=d)
        elif isinstance(xs, Union):
            for set in xs.args:
                res.extend(
                    self._map_solution(set,d)
                )
        elif isinstance(xs, Intersection):
            closure = xs.evalf()
            self._logger.warning(f'Encountered unevaluated intersection {xs} used evalf! to retrieve {closure}')
            return self._map_solution(closure,d)
        elif isinstance(xs, Reals):
            self._logger.debug(f'No unique solution.')
            res = self._purge_solutions(ds=[
                Domain(
                    start=NegativeInfinity(),
                    end=Infinity(),
                    left_open=True,
                    right_open=True
                    )
            ], d=d)
        elif isinstance(xs, SymInterval):
            self._logger.debug(f'Interval solution.')
            res = self._purge_solutions(ds=[
                Domain(
                    start=xs.start,
                    end=xs.end,
                    left_open=xs.left_open,
                    right_open=xs.right_open
                    )
            ], d=d)
        else:
            self._logger.error(xs.__class__)
            self._logger.error(xs)
            raise ValueError("Found no solution.")

        return res
    
    def _solve(self, f, d:Domain=None):
        from sympy.solvers import solveset
        from sympy import S
        from sympy.sets.sets import Interval as SymInterval
        from sympy.core.numbers import NegativeInfinity as NegInf, Infinity as Inf

        if d is None:
            boundary = S.Reals
        else:
            boundary = SymInterval(
                NegInf() if isinstance(d.start,NegativeInfinity) else ExactPymbolicToSympyMapper()(d.start),
                Inf() if isinstance(d.end,Infinity) else ExactPymbolicToSympyMapper()(d.end),
                d.left_open,
                d.right_open
            )

        self._logger.debug(f'Solving {f} within {boundary} ({d})')
        solution = solveset(f,self.symbol,boundary)
        self._logger.debug(f'Obtained {solution}, starting conversion.')

        res = self._map_solution(solution,d=d)

        return res