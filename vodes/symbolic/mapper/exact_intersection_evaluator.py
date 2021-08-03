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
from sympy import S, Interval as Domain
from sympy.series.limits import limit
from sympy.solvers import solveset

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

    def __intersection(self, f, g, b: Boundary):
        """Determine the real intersections of two functions for the given boundary. Intersections at the boundary values are ignored.

        Args:
            f: Function to intersect with g.
            g: Function to intersect with f.
            b: Boundary to limit the possible values for the intersections.

        Returns:
            __intersection: A (possibly empty) list of intersections between f and g.
        """
        res = []

        xs = solveset(
            f-g,
            domain=Domain(b.lower.value,b.upper.value,True,True)
        )

        if not xs.is_FiniteSet:
            #f=g
            return res

        for x in list(xs):
            if b.contains(x):
                # Allows equal handling of boundary and intersections
                res.append(x)
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
            intersections: The 3D array of intersections, as it was obtained using __intersections. Importantly, we ignore intersections on the boundary, as they are (per definition) equal on the boundary values and evaluations beyond that are not needed.
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
            # Closest intersection, after current intersection
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

            # Boundary was already evaluated, finished
            if res[-1][1] == min_bl:
                break

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

    def __bound_result(self, rs, b:Boundary):
        """Combines a list of bounded values to span the supplied boundary.

        Args:
            rs : List of tuples with first element being a value (e.g. expression, index) and the second element being a bounded value, describing the start of when the value is valid.

        Returns:
            __bound_result: List of 2-element tuples with first element being the value and the second element being the boundary of validity.
        """
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

        return self.__merge_bounded_results(rs=rs_bounded)
        

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
        self._logger.info(f"min: {min_res_bounded} , max: {max_res_bounded}")

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
        xs = solveset(
            f_diff,
            domain = S.Reals
        )

        if not xs.is_FiniteSet:
            self._logger.warn(f'{f_diff} was not solvable symbolically. Assuming monotonicity.')
            return res

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
            xlows = solveset(
                f_low - x,
                domain=Domain(b.lower.value,b.upper.value,True,True)
            )
            
            lower_sections = self._sboundary_sections_from_intersections(
                f= f_low,
                intersections= xlows,
                included= lambda val, extrema : val <= extrema,
                extrema=x,
                b=b
            )

            # i_up = x
            xups = solveset(
                f_up - x,
                domain=Domain(b.lower.value,b.upper.value,True,True)
            )

            upper_sections = self._sboundary_sections_from_intersections(
                f= f_up,
                intersections= xups,
                included= lambda val, extrema : val >= extrema,
                extrema=x,
                b=b
            )

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