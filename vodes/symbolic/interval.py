from functools import reduce
from abc import abstractmethod, ABC

from sys import intern

from numpy.lib.arraysetops import isin
from sympy.core.function import expand
from sympy.solvers import solveset
from vodes.symbolic.absolute import Abs
from vodes.symbolic.maximum import Max

from sympy.series.limits import limit
from vodes.symbolic.symbols import Boundary, BoundedExpression, BoundedValue, BoundedVariable
from vodes.symbolic.absolute import Abs
from vodes.symbolic.power import Power
from vodes.symbolic.bounded_mapper import BoundedMapper
from pymbolic.mapper.evaluator import EvaluationMapper

from pymbolic.mapper import RecursiveMapper
from pymbolic.polynomial import Polynomial
from pymbolic.primitives import Expression, Quotient, Variable, Sum, Product, is_constant
from pymbolic.mapper.stringifier import PREC_NONE, StringifyMapper
from pymbolic.interop.sympy import PymbolicToSympyMapper, SympyToPymbolicMapper

from sympy import solve, im, simplify, S, Eq
from typing import List
import logging

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
                elif isinstance(l, Variable) or is_constant(l):
                    lexpr = Interval(l)

                if isinstance(r, BoundedExpression):
                    bound = bound.intersect(r.bound)
                    rexpr = r.expr
                # May result from BoundedMapper
                elif isinstance(r, Variable) or is_constant(l):
                    rexpr = Interval(l)

                if not (isinstance(lexpr, Interval) and isinstance(rexpr, Interval)):
                    raise TypeError(
                        f"Cannot process {lexpr}, {rexpr} as (atleast) one entity is not an interval.")

                res.extend(op(lexpr, rexpr, bound))

        return res

    def map_constant(self, expr) -> List[Interval]:
        res = [
            BoundedExpression(
                Interval(expr),
                boundary=self._symbol.bound
            )
        ]

        self._logger.debug(f'{expr} -> {list(map(str,res))}')
        return res
  
    def map_variable(self, expr:Variable) -> List[BoundedExpression]:
        res = None
        # vexpr do not substitute the free symbol
        if not str(self._symbol.name) in self._context and self._symbol.name == expr.name:
            vexpr = Interval(expr)
        else:
            vexpr = Interval(EvaluationMapper(context=self._context)(expr))

        res = [
            BoundedExpression(
                expression=vexpr,
                boundary=self._symbol.bound
            )
        ]

        self._logger.debug(f'{expr} -> {list(map(str,res))}')
        return res

    def map_sum(self, expr:Sum) -> List[Expression]:
        res = reduce(
            lambda r, x: self.__apply(self._iadd, r, x), [
                self.rec(child) for child in expr.children
            ]
        )

        self._logger.debug(f'{expr} -> {list(map(str,res))}')
        return res

    def map_product(self, expr:Product) -> List[Expression]:
        res = reduce(
            lambda r, x: self.__apply(self._imul, r, x), [
                self.rec(child) for child in expr.children
            ]
        )

        self._logger.debug(f'{expr} -> {list(map(str,res))}')
        return res

    def map_quotient(self, expr:Quotient) -> List[Expression]:
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
        l = BoundedMapper(context=self._context, symbol=self._symbol)(expr.low)
        r = BoundedMapper(context=self._context, symbol=self._symbol)(expr.up)

        lexpr = l
        rexpr = r

        if isinstance(l, BoundedExpression):
            lexpr = l.expr

        if isinstance(r, BoundedExpression):
            rexpr = r.expr

        return [
            BoundedExpression(
                expression=Interval(lower=lexpr, upper=rexpr),
                boundary=self._symbol.bound
            )
        ]

    def map_maximum(self, expr:Max) -> List[BoundedExpression]:
        bexprs = self.rec(expr.expr)

        return [
            BoundedExpression(
                expression=item.expr.up,
                boundary=item.boundary
            ) 
            for item in bexprs
        ]

    def map_absolute(self, expr:Abs) -> List[BoundedExpression]:
        bexprs = self.rec(expr.expr)
        return [
            item for sublist in list(map(lambda bexpr : self._sabs(bexpr.expr, bexpr.bound), bexprs))
                 for item in sublist
        ]

    ##
    # Interval Arithmetic methods
    ##
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
            _idiv: A list of BoundedExpressions containing the result of the division (interval), possibly limited to subsets of the initial boundary.
        """
        rmin, rmax = r.low, r.up

        # TODO : to be handled properly
        # 0 \not in Y
        return self._imul(l, Interval(1/rmax, 1/rmin), b)
        # rmin = 0, rmax > 0
        # [1/rmax,\infty)
        # ...
        # TODO : Applied Interval Analysis

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

        # Degenerate intervals
        print(l)
        print(r)
        if lmin == lmax:
            res = Interval(lmin * rmin, lmin * rmax)
        elif rmin == rmax:
            res = Interval(lmin * rmin, lmax * rmin)

        try:
            res = Interval(
                min(lmin * rmin, lmin * rmax, lmax * rmin, lmax * rmax),
                max(lmin * rmin, lmin * rmax, lmax * rmin, lmax * rmax)
            )
        except TypeError:
            # min, max could not be handled => result of symbolic boundaries
            return self._smul(l, r, b)

        return [
            BoundedExpression(
                expression=res,
                boundary=b
            )
        ]

    ## TODO : Evaluate correctness
    def _ipow(self, l, r, b: Boundary):
        lmin, lmax = l.low, l.up
        rmin, rmax = r.low, r.up

        #if not is_constant(lmin) or not is_constant(lmax):
        return self._spow(l, r, b)

    ##
    # Symbolic Interval methods (not evalutable otherwise)
    ##
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
    def _spow(self, l:Interval, r:Interval, b: Boundary):
        pass

    @abstractmethod
    def _sabs(self, i:Interval, b: Boundary):
        pass


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


class ExactIntersectionEvaluator(ExactIntervalEvaluator):
    """Class for exactly evaluating symbolic intervals on the basis of intersections.
    TODO : Describe in depth, firstly in thesis.
    """
    def __init__(self, context: dict, symbol: Variable):
        super().__init__(context=context, symbol=symbol)

    def __intersection(self, f, g, b: Boundary):
        """Determine the real intersections of two functions for the given boundary.

        Args:
            f: Function to intersect with g.
            g: Function to intersect with f.
            b: Boundary to limit the possible values for the intersections.

        Returns:
            __intersection: A (possibly empty) list of intersections between f and g.
        """
        res = []

        xs = solveset(f-g,domain=S.Reals)

        if not xs.is_FiniteSet:
            #f=g
            return res

        for x in list(xs):
            if b.contains(x):
                # Allows equal handling of boundary and intersections
                res.append(x)
        return res

    def __intersections(self, fs:list, b: Boundary):
        """Determine the real intersections for a list of functions within the given boundary.

        Args:
            fs: Functions to evaluate.
            b: Boundary to limit the possible values for the intersections.
        
        Returns:
            __intersections: A 3D array of intersections. [i][j]=[x1,x2] -> The function i (position within fs) intersects the function j at the values x1,x2 within the boundary.
        """
        intersections = [[[] for j in range(len(fs))] for i in range(len(fs))]

        for i in range(0, len(fs)):
            for j in range(0, len(fs)):
                if i == j:
                    continue

                intersections[i][j] = self.__intersection(fs[i], fs[j], b)
                intersections[j][i] = intersections[i][j]
        
        return intersections

    def __boundary_eval(self, f, bv: BoundedValue):
        """Evaluate a function for a given bounded value.

        Args:
            f: Function to evaluate.
            bv: Value at which the function is to be evaluated. This is either done as a limit, if the value is an open bound, or exactly.

        Returns:
            __boundary_eval: The evaluation of the function f for the bounded value.
        """
        res = None

        if is_constant(f):
            return f

        # a) open -> lim
        if bv.open:
            res = limit(f, str(self._symbol), bv.value)
        # b) closed -> evaluation
        else:
            res = f.subs(str(self._symbol), bv.value)

        return res

    def __eval(self, fs: list, bv: BoundedValue, order):
        """Evaluate functions for a given bounded value and determine the extrema using the given order function.

        Args:
            fs: Functions to evaluate.
            bv: Value at which the functions are to be evaluated.
            order: Defines the comparison between function value (e.g. <, >).
        
        Returns:
            __eval: A 2-element tuple with the second element being the extreme value, as defined by the order. The first element contains the indices of the functions being evaluated to this value.
        """
        values = [None] * len(fs)

        for i in range(0, len(fs)):
            values[i] = self.__boundary_eval(fs[i], bv)

        result = ([0], values[0])
        for i in range(1, len(fs)):
            if order(values[i], result[1]):
                result = ([i], values[i])
            elif result[1] == values[i]:
                result[0].append(i)

        return result

    def __extrema_eval(self, fs: list, b: Boundary, intersections: list, order):
        """Evaluate the extrema (e.g. minimum, maximum) of the given functions within the boundary.

        Args:
            order: Determines minimum, maximum evaluation.
            fs: Functions to evaluate.
            b: Boundary within the functions are to be evaluated.
            intersections: The 3D array of intersections, as it was obtained for the functions using __intersections.

        Returns:
            __extrema_eval: The extrema function (minimum/maximum) over the boundary compared to the other functions, as index.
        """
        res = self.__eval(fs, b.lower, order)

        if len(res[0]) == 1:
            return res[0][0]

        flat_intersections = [
            item for sublist in intersections for subsublist in sublist for item in subsublist]

        bound = b.upper.value
        if len(flat_intersections) > 0:
            bound = min(flat_intersections)

        res = self.__eval(fs, BoundedValue(
            value=(b.lower.value + bound)/2,
            open=False
        ),
            order
        )

        return res[0][0]

    def __min_eval(self, fs: list, b: Boundary, intersections: list):
        """Evaluate the minimum function within the given boundary."""
        return self.__extrema_eval(fs, b, intersections, order=lambda a, b: a < b)

    def __max_eval(self, fs: list, b: Boundary, intersections: list):
        """Evaluate the maximum function within the given boundary."""
        return self.__extrema_eval(fs, b, intersections, order=lambda a, b: a > b)

    def __next_eval(self, f: int, intersections: list, b: Boundary):
        """Determine which functions are where to be evaluated next.

        Args:
            f: Current extrema function.
            intersections: The 3D array of intersections, as it was obtained using __intersections.
            b: Boundary to evaluate the next evaluation within.

        Returns:
            __next_eval: 3-element tuple. The first element is the list of functions with the next intersection with f. The second element is the intersection point with f. The third element is the closest intersection point of all functions contained within Item1.
        """
        # Item 1 : List of functions with next intersection with f
        # Item 2 : Intersection point
        # Item 3 : Closest intersection point of Item 1 with eachother
        res = [[], None, None]

        for i in range(len(intersections[f])):
            xs = [x for x in intersections[f][i] if b.contains(x)]

            if len(xs) == 0:
                continue

            if len(res[0]) == 0 or res[1] > xs[0]:
                res = [[i], xs[0], None]
            elif res[1] == xs[0]:
                res[0].append(i)
        
        if len(res[0]) == 0:
            return res

        # intersection may be tangent
        # e.g. 
        # 2.25x^{5}+0.75x^{4}-3.5x^{3}-0.5x^{2}+1.25x-0.25
        # -2.25x^{5}-5.25x^{4}-2.5x^{3}+1.5x^{2}+0.75x-0.25
        res[0].append(f)

        for i in range(len(res[0])):
            # closest intersection, after current intersection
            xs = [x for x in intersections[res[0][0]][i] if b.contains(x) and x > res[1]]

            if len(xs) == 0:
                continue
            if res[2] is None or res[2] > xs[0]:
                res[2] = xs[0]


        return res

    def __analysis(self, eval, b: Boundary, intersections: list, fs: list):
        """Evaluate the extremas (e.g. minima, maxima) of the given functions within the boundary and return a list of them with their respective boundary."""
        res = []
        extrema = None

        # (3) Determine possible switches of min/max based on intersections
        min_bl = b.lower
        max_bl = b.upper

        while True:
            candidates = []

            if extrema is None:
                # Determine the next possible intersection
                flat_intersections = [
                    item for sublist in intersections for subsublist in sublist for item in subsublist if item > min_bl.value]

                if len(flat_intersections) > 0:
                    # Evaluate up until intersection
                    max_bl = BoundedValue(
                        value=min(flat_intersections), open=True)

                candidates = [i for i in range(len(fs))]
            else:
                e_next = self.__next_eval(extrema, intersections, Boundary(
                    lower=min_bl,
                    # We start by anticipating the range to go till the end
                    upper=b.upper
                ))

                # No more intersections to handle, minimum stays minimum
                if len(e_next[0]) == 0:
                    break

                candidates = e_next[0]

                if e_next[2] is None:
                    # We evaluate until the upper boundary, as no more interesting intersection is given
                    max_bl = b.upper
                else:
                    # We evaluate up until the next intersection
                    max_bl = BoundedValue(
                        value = e_next[2],
                        open = True
                    )

            self._logger.debug(f'Evaluating {candidates} from {min_bl} to {max_bl}')
            extrema = eval([fs[i] for i in candidates], Boundary(
                lower=min_bl,
                upper=max_bl
            ), intersections)

            res.append(
                (candidates[extrema], min_bl)
            )

            self._logger.debug(f'Extrema : {fs[candidates[extrema]]} from {min_bl} to {max_bl}')

            # Start from intersection (including)
            min_bl = BoundedValue(
                value = max_bl.value,
                open = False
            )

        return res

    def __bound_result(self, rs, b:Boundary):
        rs_bounded = [None]*len(rs)

        for i in range(len(rs)):
            if (i+1) == len(rs):
                rs_bounded[i] = (rs[i][0], Boundary(
                    lower=rs[i][1],
                    upper=b.upper
                ))
            else:
                rs_bounded[i] = (rs[i][0], Boundary(
                    lower=rs[i][1],
                    upper=BoundedValue(
                        value=rs[i+1][1].value,
                        open=not rs[i+1][1].open
                    )
                ))



        rs_merged = [rs_bounded[0]]

        for i in range(len(rs_bounded) - 1):
            if rs_merged[-1][0] == rs_bounded[i+1][0]:
                rs_merged[-1][1].union(rs_bounded[i+1][1])
            else:  
                rs_merged.append(rs_bounded[i+1])

        return rs_merged

    def __sint(self, pfs:list, b:Boundary):
        """Determine the upper and lower bound of the symbolic interval and return it."""
        self._logger.debug(f"Symbolic Interval Evaluation : {list(map(str,pfs))}")
        fs = [PymbolicToSympyMapper()(f) for f in pfs]

        # (1) Determine intersection for possible min/max switch
        intersections = self.__intersections(fs,b)
        self._logger.debug(f'Intersections : {intersections}')

        min_res = self.__analysis(self.__min_eval, b, intersections, fs)
        self._logger.debug(f'Minima : {min_res}')
        max_res = self.__analysis(self.__max_eval, b, intersections, fs)
        self._logger.debug(f'Maxima : {max_res}')

        # Convert Bounded Values to Boundaries
        min_res_bounded = self.__bound_result(min_res,b)
        max_res_bounded = self.__bound_result(max_res,b)

        # TODO : Very inefficient, but hey
        result = []
        for minr in min_res_bounded:
            for maxr in max_res_bounded:
                combined = minr[1].intersect(maxr[1])

                if combined is None:
                    continue

                result.append(
                    BoundedExpression(
                        expression=Interval(
                            SympyToPymbolicMapper()(simplify(fs[minr[0]])),
                            SympyToPymbolicMapper()(simplify(fs[maxr[0]]))
                            ),
                        #Interval(pfs[minr[0]], pfs[maxr[0]]),
                        boundary=combined
                    )
                )

        self._logger.info(f"_sint: {list(map(str,pfs))}")
        self._logger.info(f"min: {min_res_bounded} , max: {max_res_bounded}")

        return result

    def _smul(self, l, r, b):
        pfs = [
            l.low * r.low,
            l.low * r.up,
            l.up * r.low,
            l.up * r.up
        ]

        return self.__sint(pfs, b)

    def __convert_to_positive(self, pf, b:Boundary):
        result = []

        f = PymbolicToSympyMapper()(pf)

        # Determine intersections with x axis
        xs = [b.lower]
        xs.extend(list(map(lambda x : BoundedValue(value=x, open=False),self.__intersection(f,0,b))))
        xs.append(b.upper)

        for i in range(len(xs) - 1):
            x = (xs[i].value + xs[i + 1].value) / 2
            y = self.__boundary_eval(f, BoundedValue(x,False))

            bi = Boundary(
                lower = xs[i],
                upper = xs[i + 1]
            )

            if y > 0:
                result.append(
                    BoundedExpression(
                        expression=pf,
                        boundary=bi
                    )
                )
            # if y == 0 -> everywhere 0 between two points
            else:
                result.append(
                    BoundedExpression(
                        expression=(-1) * pf,
                        boundary=bi
                    )
                )

        return result


    def _sabs(self, i:Interval, b:Boundary):
        lfs = self.__convert_to_positive(i.low, b)
        ufs = self.__convert_to_positive(i.up, b)

         # TODO : Very inefficient, but hey
        result = []
        for lf in lfs:
            for uf in ufs:
                combined = lf.boundary.intersect(uf.boundary)

                if combined is None:
                    continue

                    
                result.extend(self.__sint(pfs=[lf.expr,uf.expr],b=combined))

        return result


    def __get_negative(self,f,b:Boundary) -> List[Boundary]:
        """Determine where the function is negative with the given boundary, may return empty list"""

        res = []

        # x axis crossings
        # open=True -> we dont want to include x=0 portions in negative sections
        intersections = [BoundedValue(x,open=True) for x in self.__intersection(f, 0, b)]

        xs = [b.lower]
        xs.extend(intersections)
        xs.append(b.upper)

        for i in range(len(xs) - 1):
            x = (xs[i].value + xs[i+1].value)/2
            y = self.__boundary_eval(f, bv=BoundedValue(value=x,open=False))

            if y < 0:
                res.append(Boundary(
                    lower = xs[i],
                    upper = xs[i+1]
                ))
        
        return res

    def __get_positive(self,f,b:Boundary) -> List[Boundary]:
        """Determine where the function is positive with the given boundary, may return empty list"""
        return b.difference(self.__get_negative(f,b))



    def _spow_inequalities(self, pfs:list, b:Boundary):
        # The list is of the format [a^x,a^y,...]

        # TODO Try to determine min, max based on inequalities
        
        pass



    def _spow(self, l:Interval, r:Interval, b:Boundary):
        if not is_constant(r.low):
            raise ValueError("The ExactIntersectionEvaluator does not support symbolic exponents")

        if not (r.low == r.up):
            raise ValueError("The ExactIntersectionEvaluator does not support non degenerate interval exponents")

        pfs = []

        if r.low < 0:
            pfs = [
                Power(base=l.low,exponent=(-1) * r.low),
                Power(base=l.up,exponent=(-1) * r.low)
            ]
        else:
            pfs = [
                Power(base=l.low,exponent=r.low),
                Power(base=l.up,exponent=r.low)
            ]


        res = []

        # We know a (global) extrema is at l == 0
        if r.low % 2 == 0:
            lns = self.__get_negative(PymbolicToSympyMapper()(l.low),b)
            uns = self.__get_positive(PymbolicToSympyMapper()(l.up),b)

            # l.low < 0 <= l.up
            if len(lns) != 0 and len(uns) != 0:

                cbs = [ln.intersect(un) for ln in lns for un in uns]
                bns = [bn for bn in cbs if bn is not None]

                
                # none special case
                pns = b.difference(bns)

                for pn in pns:
                    res.extend(self.__sint(pfs=pfs, b=pn))

                # l.low < 0 <= l.up case
                for bn in bns:
                    # TODO : Only max to be determined
                    res.extend(self.__sint(pfs=[
                        pfs[0],
                        0,
                        pfs[1]
                    ], b=bn))
        
        if len(res) == 0:
            res = self.__sint(pfs=pfs, b=b)
        
        if r.low < 0:
            self._logger.debug(f'Encountered negative exponent, converting using division.')
            converted_res = [self._idiv(Interval(1),bexpr.expr,b=bexpr.bound) for bexpr in res]
            res = [item for sublist in converted_res for item in sublist]
            
        return res
