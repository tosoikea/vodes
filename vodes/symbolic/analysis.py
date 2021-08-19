import logging
from typing import List
from vodes.symbolic.expressions.infinity import Infinity, NegativeInfinity

# Custom Expression
from vodes.symbolic.expressions.bounded import Domain

# Custom Mapper
from vodes.symbolic.mapper.interop import ExactPymbolicToSympyMapper

# Symbolic Expression Library
from pymbolic.primitives import Expression

# Mapper
from pymbolic.interop.sympy import SympyToPymbolicMapper

class Analysis:
    def __init__(self, pf:Expression) -> None:
        assert not (pf is None)
        self.__pf = pf
        self.__f = ExactPymbolicToSympyMapper()(self.__pf)

        if len(self.__f.free_symbols) > 1:
            raise ValueError("Not supporting functions with more than one free variable. Consider encoding the parameter as a vector.")

        if len(self.__f.free_symbols) == 1:
            self.__symbol = list(self.__f.free_symbols)[0]
        else:
            self.__symbol = None

        self._logger = logging.getLogger(__name__)


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

    def equals(self, y, d:Domain) -> List[Domain]:
        ask = self.__f - y
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

    def taylor(self, a, n:int=1):
        from sympy import series
        assert n > 0

        if isinstance(a,Expression):
            a = ExactPymbolicToSympyMapper()(a)

        # (1) Construct the taylor polynomial
        f = series(self.func,self.symbol,x0=a,n=n)

        bigO = f.getO()

        f = f.removeO()

        # (2) Use the taylor remainder estimation theorem
        # TODO

        return SympyToPymbolicMapper()(f)

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
        from sympy.core.basic import Basic
        from sympy import CRootOf
        
        if isinstance(d.start,CRootOf):
            raise ValueError(f'Encountered unsolved root {d.start}')
        
        if isinstance(d.end,CRootOf):
            raise ValueError(f'Encountered unsolved root {d.start}')

        if isinstance(d.start,Basic):
            d.start = SympyToPymbolicMapper()(d.start)
        
        if isinstance(d.end,Basic):
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
            self._logger.debug(_d)
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
        from sympy.sets.sets import Union, Intersection
        from sympy.sets.fancysets import Reals
        from sympy import S

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
        else:
            self._logger.error(xs.__class__)
            self._logger.error(xs)
            raise ValueError("Found no solution.")

        return res
    
    def _solve(self, f, d:Domain=None):
        from sympy.solvers import solveset
        from sympy import S
    
        res = self._map_solution(xs=solveset(f,domain=S.Reals),d=d)
        self._logger.debug(f'{f} -> {res}')

        return res