from functools import reduce
from abc import abstractmethod, ABC

from sys import intern
from vodes.symbolic.absolute import AbsoluteEvaluator
from sympy.geometry.util import intersection
from sympy.matrices.expressions.matadd import combine

from sympy.series.limits import limit
from vodes.symbolic.symbols import Boundary, BoundedExpression, BoundedValue, BoundedVariable
from vodes.symbolic.absolute import Abs

from pymbolic.mapper import RecursiveMapper
from pymbolic.polynomial import Polynomial
from pymbolic.primitives import Expression, Lookup, Quotient, Variable, Sum, Product
from pymbolic.mapper.evaluator import EvaluationMapper
from pymbolic.mapper.stringifier import PREC_NONE, StringifyMapper
from pymbolic.interop.sympy import PymbolicToSympyMapper, SympyToPymbolicMapper

from sympy import solve, im
from typing import List

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
        assert(lower)

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

    def map_approx_polynomial(self, expr, enclosing_prec, *args, **kwargs):
        return str(Polynomial(expr.Base, expr.Data))


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
        assert(context)
        assert(symbol)

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
                    bound.intersect(l.bound)
                    lexpr = l.expr

                if isinstance(r, BoundedExpression):
                    bound.intersect(r.bound)
                    rexpr = r.expr

                if not (isinstance(lexpr, Interval) and isinstance(rexpr, Interval)):
                    raise TypeError(
                        f"Cannot process {lexpr}, {rexpr} as (atleast) one entity is not an interval.")

                res.extend(op(lexpr, rexpr, bound))

        return res

    def map_constant(self, expr) -> List[Interval]:
        return [
            Interval(expr)
        ]

    def map_variable(self, expr:Variable) -> List[BoundedExpression]:
        res = None
        # we do not substitute the free symbol
        if not str(self._symbol.expr) in self._context and self._symbol == expr:
            res = Interval(expr)
        else:
            res = EvaluationMapper(context=self._context)(expr)

        return [
            BoundedExpression(
                Interval(res),
                boundary=self._symbol.bound
            )
        ]

    def map_sum(self, expr:Sum) -> List[Expression]:
        return reduce(
            lambda r, x: self.__apply(self._iadd, r, x), [
                self.rec(child) for child in expr.children
            ]
        )

    def map_product(self, expr:Product) -> List[Expression]:
        return reduce(
            lambda r, x: self.__apply(self._imul, r, x), [
                self.rec(child) for child in expr.children
            ]
        )

    def map_quotient(self, expr:Quotient) -> List[Expression]:
        return self.__apply(
            self._idiv,
            self.rec(expr.numerator),
            self.rec(expr.denominator)
        )

    def map_interval(self, expr) -> List[BoundedExpression]:
        return [
            BoundedExpression(
                expression=expr,
                boundary=self._symbol.bound
            )
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

        # TODO : 0 \in Y assumed -> to be handled
        return self._imul(l, Interval(1/rmin, 1/rmax), b)

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

    #
    # TODO : Handle special cases (e.g. 0**-1)
    #
    def _ipow(self, l, r, b: Boundary):
        res = None
        lmin, lmax = l.low, l.up
        rmin, rmax = r.low, r.up

        try:
            res = Interval(
                min(lmin ** rmin, lmin ** rmax, lmax ** rmin, lmax ** rmax),
                max(lmin ** rmin, lmin ** rmax, lmax ** rmin, lmax ** rmax)
            )
        except TypeError:
            return self._spow(l, r, b)

        return [
            BoundedExpression(
                expression=res,
                boundary=b
            )
        ]

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

        for x in solve(f - g):
            if im(x) == 0 and b.contains(x):
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
        intersections = [[[]] * len(fs)] * len(fs)

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
            return res

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

    def __next_eval(self, f: int, intersections: list, b: int):
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
            xs = [x for x in intersections[f][i] if x > b]

            if len(xs) == 0:
                continue

            if len(res[0]) == 0 or res[1] > xs[0]:
                res = ([i], xs[0])
            elif res[1] == xs[0]:
                res[0].append(i)

        for i in range(len(res[0])):
            xs = [x for x in intersections[res[0][0]][i] if x > b]

            if len(xs) == 0:
                continue

            if res[2] is None or res[2] > xs[0]:
                res[2] = xs[0]

        return res

    def __analysis(self, eval, b: Boundary, intersections: list, fs: list):
        """Evaluate the extremas (e.g. minima, maxima) of the given functions within the boundary and return a list of them with their respective boundary."""
        res = []
        extrema = None
        min_b = None

        # (3) Determine possible switches of min/max based on intersections
        while True:
            min_bl = b.lower
            max_bl = b.upper
            candidates = []

            if extrema is None:
                flat_intersections = [
                    item for sublist in intersections for subsublist in sublist for item in subsublist]

                if len(flat_intersections) > 0:
                    max_bl = BoundedValue(
                        value=min(flat_intersections), open=False)

                candidates = fs
            else:
                e_next = self.__next_eval(extrema, intersections, min_b)

                # No more intersections to handle, minimum stays minimum
                if len(e_next[0]) == 0:
                    break

                candidates = e_next[0]

            extrema = eval(candidates, Boundary(
                lower=min_bl,
                upper=max_bl
            ), intersections)

            res.append(
                (extrema, min_bl)
            )

            min_b = min_bl.value

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

        return rs_bounded

    def __sint(self, pfs:list, b:Boundary):
        """Determine the upper and lower bound of the symbolic interval and return it."""
        fs = [PymbolicToSympyMapper()(f) for f in pfs]
        print(fs)

        # (1) Determine intersection for possible min/max switch
        intersections = self.__intersections(fs,b)

        min_res = self.__analysis(self.__min_eval, b, intersections, fs)
        max_res = self.__analysis(self.__max_eval, b, intersections, fs)

        # Convert Bounded Values to Boundaries
        min_res_bounded = self.__bound_result(min_res,b)
        max_res_bounded = self.__bound_result(max_res,b)

        print(min_res_bounded)
        print(max_res_bounded)

        # TODO : Very inefficient, but hey
        result = []
        for minr in min_res_bounded:
            for maxr in max_res_bounded:
                combined = minr[1].intersect(maxr[1])

                if combined is None:
                    continue

                result.append(
                    BoundedExpression(
                        expression=Interval(pfs[minr[0]], pfs[maxr[0]]),
                        boundary=combined
                    )
                )

        return result

    def _smul(self, l, r, b):
        pfs = [
            l.low * r.low,
            l.low * r.up,
            l.up * r.low,
            l.up * r.up
        ]

        return self.__sint(pfs, b)



    def _sabs(self, i:Interval, b:Boundary):
        pfs = [
            AbsoluteEvaluator(absolute=True,b=b)(i.up),
            AbsoluteEvaluator(absolute=True,b=b)(i.low)
        ]

        return self.__sint(pfs, b)

    def _spow(self, l, r):
        pass
