from functools import reduce
from abc import abstractmethod, ABC

from sys import intern
from sympy.geometry.util import intersection
from sympy.matrices.expressions.matadd import combine

from sympy.series.limits import limit
from vodes.symbolic.symbols import Boundary, BoundedExpression, BoundedValue, BoundedVariable

from pymbolic.mapper import RecursiveMapper
from pymbolic.polynomial import Polynomial
from pymbolic.primitives import Expression, Lookup, Variable
from pymbolic.mapper.evaluator import EvaluationMapper
from pymbolic.mapper.stringifier import PREC_NONE, StringifyMapper
from pymbolic.interop.sympy import PymbolicToSympyMapper, SympyToPymbolicMapper

from sympy import solve, im
##
# Auxiliary interval functions.
# src : https://link.springer.com/content/pdf/10.1007/3-540-36599-0_7.pdf
# src : https://www.math.kit.edu/ianm2/~kulisch/media/arjpkx.pdf
###


class Interval(Expression):
    init_arg_names = ("lower", "upper",)

    def __init__(self, lower, upper=None):
        # degenerate interval
        if upper is None:
            upper = lower

        self.lower = lower
        self.upper = upper

    def __getinitargs__(self):
        return self.lower, self.upper

    @property
    def low(self):
        return self.lower

    @property
    def up(self):
        return self.upper

    def make_stringifier(self, originating_stringifier=None):
        return IntervalStringifyMapper()

    mapper_method = intern("map_interval")

    # Inequality
    def __abs__(self):
        return max(abs(self.low), abs(self.up))


class IntervalStringifyMapper(StringifyMapper):
    def map_interval(self, expr, enclosing_prec, *args, **kwargs):
        lower = self.rec(expr.low, PREC_NONE, *args, **kwargs)
        upper = self.rec(expr.up, PREC_NONE, *args, **kwargs)

        return f'[{lower},{upper}]'

    def map_approx_polynomial(self, expr, enclosing_prec, *args, **kwargs):
        return str(Polynomial(expr.Base, expr.Data))


class SymbolicIntervalEvaluator(ABC, RecursiveMapper):
    def __init__(self, context: dict, symbol: BoundedVariable):
        self.context = context
        self.symbol = symbol

    def __apply(self, op, ls: list, rs: list):
        res = []

        for l in ls:
            for r in rs:
                bound = self.symbol.bound
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

    def map_constant(self, expr):
        return [
            Interval(expr)
        ]

    def map_variable(self, expr):
        res = None
        # we do not substitute the free symbol
        if not str(self.symbol.expr) in self.context and self.symbol == expr:
            res = Interval(expr)
        else:
            res = EvaluationMapper(context=self.context)(expr)

        return [
            BoundedExpression(
                Interval(res),
                boundary=self.symbol.bound
            )
        ]

    def map_sum(self, expr):
        return reduce(
            lambda r, x: self.__apply(self._iadd, r, x), [
                self.rec(child) for child in expr.children
            ]
        )

    def map_product(self, expr):
        return reduce(
            lambda r, x: self.__apply(self._imul, r, x), [
                self.rec(child) for child in expr.children
            ]
        )

    def map_quotient(self, expr):
        return self.__apply(
            self._idiv,
            self.rec(expr.numerator),
            self.rec(expr.denominator)
        )

    def map_interval(self, expr):
        return [
            BoundedExpression(
                expression=expr,
                boundary=self.symbol.bound
            )
        ]

    # TODO : Only fabs supported
    def map_call(self, expr, *args):
        def make_f(name):
            return Lookup(Variable("math"), name)

        if expr.function == make_f('fabs'):
            children = [self.rec(par, *args)
                        for i, par in enumerate(expr.parameters)]
            return abs(*children)
        else:
            raise NotImplementedError(
                f"The function {expr.function} is not supported")

    ##
    # Interval Arithmetic methods
    ##
    def _iadd(self, l, r, b: Boundary):
        return [
            BoundedExpression(
                expression=Interval(l.low + r.low, l.up + r.up),
                boundary=b
            )
        ]

    def _isub(self, l, r, b: Boundary):
        return [
            BoundedExpression(
                expression=Interval(l.low - r.up, l.up - r.low),
                boundary=b
            )
        ]

    def _idiv(self, l, r, b: Boundary):
        rmin, rmax = r.lower, r.upper

        # TODO : 0 \in Y assumed -> to be handled
        return self._imul(l, Interval(1/rmin, 1/rmax), b)

    def _imul(self, l, r, b: Boundary):
        res = None
        lmin, lmax = l.lower, l.upper
        rmin, rmax = r.lower, r.upper

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
        lmin, lmax = l.lower, l.upper
        rmin, rmax = r.lower, r.upper

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
    def _smul(self, l, r, b: Boundary):
        pass

    @abstractmethod
    def _spow(self, l, r, b: Boundary):
        pass


class ApproxIntervalEvaluator(SymbolicIntervalEvaluator):
    def __init__(self, context: dict, symbol: Variable):
        super().__init__(context=context, symbol=symbol)


class ExactIntervalEvaluator(SymbolicIntervalEvaluator):
    def __init__(self, context: dict, symbol: Variable):
        super().__init__(context=context, symbol=symbol)


class ExactIntersectionEvaluator(ExactIntervalEvaluator):
    def __init__(self, context: dict, symbol: Variable):
        super().__init__(context=context, symbol=symbol)

    def __intersection(self, f, g, b: Boundary):
        res = []

        for x in solve(f - g):
            if im(x) == 0 and b.contains(x):
                # Allows equal handling of boundary and intersections
                res.append(x)

        return res

    def __boundary_eval(self, f, bv: BoundedValue):
        res = None

        # a) open -> lim
        if bv.open:
            res = limit(f, str(self.symbol), bv.value)
        # b) closed -> evaluation
        else:
            res = f.subs(str(self.symbol), bv.value)

        return res

    def __eval(self, fs: list, b: BoundedValue, order):
        values = [None] * len(fs)

        for i in range(0, len(fs)):
            values[i] = self.__boundary_eval(fs[i], b)

        result = ([0], values[0])
        for i in range(1, len(fs)):
            if order(values[i], result[1]):
                result = ([i], values[i])
            elif result[1] == values[i]:
                result[0].append(i)

        return result

    def __extrema_eval(self, fs: list, b: Boundary, intersections: list, order):
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

        return res

    def __min_eval(self, fs: list, b: Boundary, intersections: list):
        return self.__extrema_eval(fs, b, intersections, order=lambda a, b: a < b)

    def __max_eval(self, fs: list, b: Boundary, intersections: list):
        return self.__extrema_eval(fs, b, intersections, order=lambda a, b: a > b)

    def __next_eval(self, f: int, intersections: list, b: int):
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
            ), intersections)[0][0]

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


    def _smul(self, l, r, b):
        pfs = [
            l.low * r.low,
            l.low * r.up,
            l.up * r.low,
            l.up * r.up
        ]

        fs = [PymbolicToSympyMapper()(f) for f in pfs]
        print(fs)


        # (1) Determine intersection for possible min/max switch
        intersections = [[[]] * len(fs)] * len(fs)

        for i in range(0, len(fs)):
            for j in range(0, len(fs)):
                if i == j:
                    continue

                intersections[i][j] = self.__intersection(fs[i], fs[j], b)
                intersections[j][i] = intersections[i][j]

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

    def _spow(self, l, r):
        pass
