from typing import List
from vodes.symbolic import symbols

from vodes.symbolic.interval import Interval, ExactIntervalEvaluator
from vodes.symbolic.symbols import Boundary, BoundedExpression, BoundedValue
from vodes.symbolic.mapper.interop import ExactPymbolicToSympyMapper

# Expression library
from pymbolic.primitives import Variable, Power, Expression
from pymbolic.mapper.differentiator import DifferentiationMapper as DM
from pymbolic.interop.sympy import SympyToPymbolicMapper

# Symbolic library for intersections, solving...
from sympy.series.limits import limit

class DifferentiationMapper(DM):
    def map_variable(self, expr:Variable, *args):
        if expr.name == self.variable.name:
            return 1
        else:
            return 0

class ExactIntersectionEvaluator(ExactIntervalEvaluator):
    """Class for exactly evaluating symbolic intervals on the basis of intersections.
    TODO : Describe in depth, firstly in thesis.
    """
    def __init__(self, context: dict, symbol: Variable):
        super().__init__(context=context, symbol=symbol)

    def __solve(self, f, b:Boundary=None):
        from sympy.solvers import solveset
        from sympy import im
        res= []
        xs = solveset(f)

        if not xs.is_FiniteSet:
            self._logger.info(f'{f} has no unique solution. Handling, as if no solution present.')
            return res

        for x in xs:
            #TODO : Make real/complex support a parameter choice
            if im(x) != 0:
                self._logger.debug(f'Dropping {x}.')
                continue

            if (b is None) or b.contains(x):
                res.append(x)

        return res

    def __intersection(self, f, g, b: Boundary):
        """Determine the real intersections of two functions for the given boundary. Intersections at the boundary values are ignored.

        Args:
            f: Function to intersect with g.
            g: Function to intersect with f.
            b: Boundary to limit the possible values for the intersections.

        Returns:
            __intersection: A (possibly empty) list of intersections between f and g.
        """
        res = self.__solve(
            f=f-g,
            b=Boundary(
                lower=BoundedValue(value=b.lower.value, open=True),
                upper=BoundedValue(value=b.upper.value, open=True)
            )
        )

        print(res)
    
        return res

    def __intersections(self, fs:list, b: Boundary):
        """Determine the real intersections for a list of functions within the given boundary. Intersections at the boundary values are ignored.

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

        # is_costant() of sympy objects
        if f.is_constant():
            return f

        # a) open -> lim
        if bv.open:
            res = limit(f, str(self._symbol), bv.value)
        # b) closed -> evaluationDifferentiationMapper
        else:
            res = f.subs(str(self._symbol), bv.value)

        return res

    def __min_eval(self, x0, x1,fs:list):
        """Evaluate the minimum function within the given boundary."""
        return self.__eval(fs=fs,x0=x0, x1=x1, order=lambda a, b: a < b)

    def __max_eval(self, x0, x1,fs:list):
        """Evaluate the maximum function within the given boundary."""
        return self.__eval(fs=fs,x0=x0, x1=x1, order=lambda a, b: a > b)

    def __eval(self, order, x0, x1, fs:list):
        """Evaluate functions between two values and determine the extrema using the given order function.

        Args:
            fs: Functions to evaluate.
            x0: Start of the section
            x1: End of the section
            order: Defines the comparison between function value (e.g. <, >).
        
        Returns:
            __eval: 
        """

        values = [None] * len(fs)
        x = x0 + (x1-x0)/2

        self._logger.debug(f'Evaluating {fs} at {x} ({x0}, {x1})')

        for i in range(0, len(fs)):
            values[i] = self.__boundary_eval(fs[i], BoundedValue(
                value=x,
                open=False
            ))
        
        result = (0, values[0])
        for i in range(1, len(fs)):
            if order(values[i], result[1]):
                result = (i, values[i])
            elif result[1] == values[i]:
                # As the evaluation takes places between two intersections, the functions must be equivalent, as only these cases of identical values are not contained within the intersections.
                self._logger.info(f'Encountered equivalent functions {fs[result[0]]} and {fs[i]} within the section.')

        return result[0]


    def __analysis(self, eval, b: Boundary, intersections: list, fs: list):
        """Evaluate the extremas (e.g. minima, maxima) of the given functions within the boundary and return a list of them with their respective boundary."""
        res = []

        lower = b.lower
        candidates = [i for i in range(len(fs))]
        
        # There is still domain left to analyze
        while lower.value < b.upper.value:
            x0 = lower.value

            # Determine the next possible intersection
            xs = [item for sublist in intersections for subsublist in sublist for item in subsublist if item > x0]
            xs.append(b.upper.value)

            x1 = min(xs)

            # Determine extremum
            extr_i = candidates[
                eval(
                    x0=x0,
                    x1=x1,
                    fs=[fs[i] for i in candidates]
                    )
                ]

            # Determine sub-domain
            upper = b.upper
            for i in range(len(intersections[extr_i])):
                xs = [item for item in intersections[extr_i][i] if item > x0 and item <= upper.value]

                if len(xs) == 0:
                    continue
                else:
                    t = min(xs)

                    if t < upper.value:
                        # may be tangent
                        candidates = [extr_i, i]

                        upper = BoundedValue(value=t,open=True)
                    else:
                        candidates.append(i)
                    
                    
            self._logger.debug(f'Extrema : {fs[extr_i]} from {lower} to {upper}')

            res.append(
                (extr_i, Boundary(lower=lower,upper=upper))
            )

            lower = BoundedValue(value=upper.value, open=not upper.open)

        return res


    def __merge_bounded_results(self, rs):
        """Merges bounds based on equivalency of bounded value.
        
        Args:
            rs : List of tuples with first element being the expression (index, identity value ...) and the second element being the boundary.

        Returns:
            __merge_bounded_results: List of 2-element tuples."""

        rs_merged = [rs[0]]

        for i in range(len(rs) - 1):
            if rs_merged[-1][0] == rs[i+1][0]:
                rs_merged[-1] = (
                    rs_merged[-1][0],
                    rs_merged[-1][1].union(rs[i+1][1])
                )
            else:  
                rs_merged.append(rs[i+1])

        return rs_merged

    def __merge_bounded_expressions(self, bexprs):
        """Merges bounds based on equivalency of expression"""
        return list(
            map(
                lambda tple: BoundedExpression(expression=tple[0],boundary=tple[1]),
                self.__merge_bounded_results(list(map(lambda bexpr : (bexpr.expr,bexpr.bound),bexprs)))
            ) 
        )

    def __sint(self, pfs:list, b:Boundary):
        """Determine the upper and lower bound of the symbolic interval and return it."""
        self._logger.debug(f"Symbolic Interval Evaluation : {list(map(str,pfs))}")

        fs = [ExactPymbolicToSympyMapper()(f) for f in pfs]

        # (1) Determine intersection for possible min/max switch
        intersections = self.__intersections(fs,b)
        self._logger.debug(f'Intersections : {intersections}')

        min_res = self.__analysis(self.__min_eval, b, intersections, fs)
        self._logger.debug(f'Minima : {min_res}')
        max_res = self.__analysis(self.__max_eval, b, intersections, fs)
        self._logger.debug(f'Maxima : {max_res}')

        # Convert Bounded Values to Boundaries

        # TODO : Very inefficient, but hey
        result = []
        for minr in min_res:
            for maxr in max_res:
                combined = minr[1].intersect(maxr[1])

                if combined is None:
                    continue

                result.append(
                    BoundedExpression(
                        expression=Interval(
                            pfs[minr[0]],
                            pfs[maxr[0]]
                            ),
                        # Simplifies the expressions, however quotients are converted to floating numbers!
                        #Interval(
                        #   SympyToPymbolicMapper()(simplify(fs[minr[0]])),
                        #   SympyToPymbolicMapper()(simplify(fs[maxr[0]]))
                        #),
                        boundary=combined
                    )
                )

        self._logger.info(f"_sint: {list(map(str,pfs))}")
        self._logger.info(f"min: {min_res} , max: {max_res}")

        return result

    def __convert_to_positive(self, pf:Expression, b:Boundary):
        """Determines the absolute expression.
        
        Args:
            pf (Expression): The expression to convert to an expression allowing for only positive (absolute) values.
            b (Boundary): Boundary within the expression is valid and to be converted.
            
        Returns:
            __convert_to_positive: A list of BoundedExpressions, which in sum describe the absolute expression."""
        result = []

        f = ExactPymbolicToSympyMapper()(pf)

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

    ### POWER
    def _sboundary_sections_from_intersections(self, f, intersections:list, included, extrema, b:Boundary):
        """Convert the intersections of the boundary expression with an extrema value into a sectioning within the boundary. Importantly, we ignore intersections on the boundary, as they are (per definition) equal on the boundary values and evaluations beyond that are not needed."""
        res = [
            # (
            #   Boundary ( start, end ),
            #   [ extrema ]
            # )
        ]

        lower = b.lower

        # 1. Determine sectioning
        for x in intersections:
            if not b.contains(x):
                continue
            
            bound = Boundary(
                lower = lower,
                upper = BoundedValue(value=x, open=False if x != b.upper.value else b.upper.open)
            )

            res.append(
                (bound, [])
            )

            lower = BoundedValue(value=bound.upper.value, open=not bound.upper.open)

        if len(res) == 0 or res[-1][0].upper.value < b.upper.value:
            bound = Boundary(
                lower = lower,
                upper = b.upper
            )

            res.append(
                (bound, [])
            )

        # 2. Determine if extrema is given within sections
        for section in res:
            val = self.__boundary_eval(
                f, 
                BoundedValue(
                    value=(section[0].lower.value + section[0].upper.value)/2,
                    open=False
                )
            )
            
            is_included = included(val, extrema)

            if is_included:
                section[1].append(SympyToPymbolicMapper()(extrema))

        return res

    def _sboundary_sections(self, f_diff, f_low, f_up, b:Boundary):
        from sympy import im
        res = [
            # (
            #   Boundary ( start, end ),
            #   [ extrema ]
            # )
            (
                b,
                []
            )
        ]

        # 1. Determine zero points
        # TODO : Evaluate constraint, if needed
        xs = self.__solve(f = f_diff)

        # TODO : e.g. x^3 -> differentiation between global/local extrema required
        # TODO : Ensure sorting

        # The values within the set xs ensure, that f'(x) = 0 and indicate a (local) extrema.
        # We now search our boundaries for intersections with this value.
        # If no intersection is given, the value is always within the boundary.
        # Otherwise, the value is only within certain value ranges encluded in the interval.


        # 2. Determine sectioning for each zero point
        for x in xs:
            self._logger.debug(f'Analyzing zero point {x}.')

            # i_low = x
            xlows = self.__solve(
                f=f_low-x,
                b=Boundary(
                    lower=BoundedValue(value=b.lower.value, open=True),
                    upper=BoundedValue(value=b.upper.value, open=True)
                )
            )
            
            lower_sections = self._sboundary_sections_from_intersections(
                f= f_low,
                intersections= xlows,
                included= lambda val, extrema : val <= extrema,
                extrema=x,
                b=b
            )

            # i_up = x
            xups = self.__solve(
                f=f_up - x,
                b=Boundary(
                    lower=BoundedValue(value=b.lower.value, open=True),
                    upper=BoundedValue(value=b.upper.value, open=True)
                )
            )

            upper_sections = self._sboundary_sections_from_intersections(
                f= f_up,
                intersections= xups,
                included= lambda val, extrema : val >= extrema,
                extrema=x,
                b=b
            )

            print(xlows)
            print(lower_sections)
            print(xups)
            print(upper_sections)

            # Merge with current res
            t_res = []

            for lower_section in lower_sections:
                for upper_section in upper_sections:
                    combined = lower_section[0].intersect(upper_section[0])

                    if combined is None:
                        continue

                    included = lower_section[1]

                    # One includes, other excludes => In summary excluded
                    if len(lower_section[1]) != len(upper_section[1]):
                        included = []

                    for existing_section in res:
                        existing_combined = combined.intersect(existing_section[0])

                        if existing_combined is None:
                            continue

                        included.extend(existing_section[1])

                        t_res.append(
                            (
                                existing_combined,
                                included
                            )
                        )

            res = t_res

        #TODO : Minimize
        self._logger.debug(f'Determined the following extrema inclusions. {res}')
        return res



    ####
    # SYMBOLIC INTERVAL INTERFACE
    ####
    def _sabs(self, i:Interval, b:Boundary):
        """Absolute function in the case of symbolic interval boundaries. This class uses intersections to determine the exact min/max expressions.
        
        Args:
            i (Interval): The interval to apply the absolute function on.
            b (Boundary): The boundary ('Gültigkeitsbereich') of the absolute function application.

        Returns:
            _sabs: A list of BoundedExpressions containing the result of applying the absolute function on the symbolic interval, possibly limited to subsets of the initial boundary.
        """
        lfs = self.__convert_to_positive(i.low, b)
        ufs = self.__convert_to_positive(i.up, b)

         # TODO : Improve efficiency
        result = []
        for lf in lfs:
            for uf in ufs:
                combined = lf.boundary.intersect(uf.boundary)

                if combined is None:
                    continue

                    
                result.extend(self.__sint(pfs=[lf.expr,uf.expr],b=combined))

        return result
        
    def _smul(self, l, r, b):
        """Interval multiplication in the case of symbolic interval boundaries. This class uses intersections to determine the exact min/max expressions.

        Args:
            l (Interval): The left parameter of the multiplication.
            r (Interval): The right parameter of the multiplication.
            b (Boundary): The boundary ('Gültigkeitsbereich') of the multiplication.

        Returns:
            _smul: A list of BoundedExpressions containing the result of the symbolic multiplication (interval), possibly limited to subsets of the initial boundary.
        """
        pfs = [
            l.low * r.low,
            l.low * r.up,
            l.up * r.low,
            l.up * r.up
        ]

        return self.__sint(pfs, b)
    
    def _spow(self, l:Interval, r:Interval, b:Boundary):
        if not self.is_constant(r.low):
            raise ValueError("The ExactIntersectionEvaluator does not support symbolic exponents")

        if not (r.low == r.up):
            raise ValueError("The ExactIntersectionEvaluator does not support non degenerate interval exponents")

        res = []
        
        # 1. Differentiate the Operator
        f = lambda x : Power(x,r.low)
        f_diff = DifferentiationMapper(Variable("x"))(
            Power(
                Variable("x"),
                r.low
            )
        )

        self._logger.debug(f'Using {f_diff} for determining the extrema of the power operator.')

        # 2. Determine sections
        sym_expr = ExactPymbolicToSympyMapper()(f_diff)
        sections = self._sboundary_sections(
            f_diff = sym_expr,
            f_low = ExactPymbolicToSympyMapper()(l.low),
            f_up =  ExactPymbolicToSympyMapper()(l.up),
            b=b
        )

        # 3. Determine min, max for sections
        exprs = [
            f(l.low),
            f(l.up)
        ]

        # TODO Determine monoton sections
        # Function can then simply be applied to boundaries
        # https://epubs.siam.org/doi/pdf/10.1137/1.9780898717716.ch5
        for section in sections:
            pfs = set(exprs)

            # All extrema can be possible boundaries
            # e.g. [-x,x] ** 2 
            # != [(-x)**2,(x)**2]
            # == [min(0,(-x)**2,(x)**2),max(0,(-x)**2,(x)**2)] = [0,x**2]
            for extrema in section[1]:
                # TODO : Here, we potentially add floating numbers, this makes my face turn this way :(
                # One could round to a "near" rational number. However, there is always going to be a "closer" rational number between them.
                #    for x in xs:
                #        if not b.contains(x):
                #            continue
                #            
                #        pfs.add(evaluate(substitute(expr,{self._symbol.name: SympyToPymbolicMapper()(x)})))
                pfs.add(f(extrema))

            self._logger.debug(f'Determining min/max for {pfs}.')
            
            res.extend(
                self.__sint(
                    pfs=list(pfs),
                    b=section[0]
                )
            )
       
        return self.__merge_bounded_expressions(res)