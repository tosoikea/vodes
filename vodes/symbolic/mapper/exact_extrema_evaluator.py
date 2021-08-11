
from abc import abstractmethod, ABC

from vodes.symbolic.interval import Interval, ExactIntervalEvaluator
from vodes.symbolic.symbols import Boundary, BoundedValue, BoundedExpression
from vodes.symbolic.power import Power

# Expression library
from pymbolic.primitives import Variable, Expression

class ExactExtremaEvaluator(ExactIntervalEvaluator, ABC):
    """Class for determining the exact boundaries of intervals on the basis of function analysis."""

    def __init__(self, context: dict, symbol: Variable):
        super().__init__(context=context, symbol=symbol)

    @abstractmethod
    def _get_symbolic_interval(self, pfs:list, b: Boundary) -> list:
        """Construct symbolic intervals based on the included expressions within the given boundary. Possibly returns a list of multiple intervals split into sub-boundaries."""
        pass

    def symbolic_expression(self, expr:Expression, iv:Interval, extrema:list, b:Boundary) -> list:
        """
        TODO : Documentation
        """
        from pymbolic import substitute
        res = []

        # 1. Determine Sections
        sections = self.__get_sectioning(iv=iv, extrema=extrema,b=b)

        # 2. Determine min, max for sections
        exprs = [
            substitute(expr,{Variable("x"): iv.low}),
            substitute(expr,{Variable("x"): iv.up})
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
                pfs.add(substitute(expr,{Variable("x"): extrema}))

            self._logger.debug(f'Determining min/max for {pfs}.')

            res.extend(
                self._get_symbolic_interval(
                    pfs=list(pfs),
                    b=section[0]
                )
            )

        print(res)
        return self.__merge_bounded_expressions(res)

    def __get_sectioning(self,iv:Interval,extrema:list,b:Boundary) -> list:
        """Function to determine the sectioning based on when the interval contains the listed extremas."""
        from vodes.symbolic.mapper.interop import ExactPymbolicToSympyMapper

        # The values within the extremas ensure that f'(x) = 0 and indicate a (local) extrema.
        # We now search our boundaries for intersections with this value.
        # If no intersection is given, the value is always within the boundary.
        # Otherwise, the value is only within certain value ranges encluded in the interval.
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

        f_low = ExactPymbolicToSympyMapper()(iv.low)
        f_up =  ExactPymbolicToSympyMapper()(iv.up)

        for x in extrema:
            self._logger.debug(f'Analyzing zero point {x}.')
            xb = Boundary(
                    lower=BoundedValue(value=b.lower.value, open=True),
                    upper=BoundedValue(value=b.upper.value, open=True)
                )

            # i_low = x
            xlows = self._solve(
                f=f_low-x,
                b=xb
            )
            
            lower_sections = self.__get_boundary_sections(
                f= f_low,
                intersections= xlows,
                included= lambda val, extrema : val <= extrema,
                extrema=x,
                b=b
            )

            # i_up = x
            xups = self._solve(
                f=f_up - x,
                b=xb
            )

            upper_sections = self.__get_boundary_sections(
                f= f_up,
                intersections= xups,
                included= lambda val, extrema : val >= extrema,
                extrema=x,
                b=b
            )

            res = self.__merge_sections(
                existing=res,
                expansions=self.__get_interval_sections(
                    upper_sections=upper_sections,
                    lower_sections=lower_sections,
                    b=b
                )
            )

        #TODO : Minimize
        self._logger.debug(f'Determined the following extrema inclusions. {res}')
        return res

    def __get_boundary_sections(self, f, intersections:list, included, extrema, b:Boundary):
        """Convert the intersections of the boundary expression with an extrema value into a sectioning within the boundary. Importantly, we ignore intersections on the boundary, as they are (per definition) equal on the boundary values and evaluations beyond that are not needed."""
        from pymbolic.interop.sympy import SympyToPymbolicMapper

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
            val = self._boundary_eval(
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

    def __get_interval_sections(self,upper_sections:list,lower_sections:list,b:Boundary):
        """Determine the inclusion of an extrema (section) based on the inclusion of the upper and lower boundaries. For the interval to include an extrema, both boundaries have to include it."""        
        res = []

        for lower_section in lower_sections:
            for upper_section in upper_sections:
                combined = lower_section[0].intersect(upper_section[0])

                if combined is None:
                    continue

                included = lower_section[1]

                # One includes, other excludes => In summary excluded
                # May only include one extrema, as only one extrema was tested for inclusion
                if len(lower_section[1]) != len(upper_section[1]):
                    included = []

                res.append(
                    (
                        combined,
                        included
                    )
                )
        self._logger.debug(f'Combined boundary section {lower_sections}(l) and {upper_sections}(u) to {res}')
        return res

    ####
    # SYMBOLIC INTERVAL INTERFACE
    ####
    def _spow(self, l:Interval, r:Interval, b:Boundary):
        """Power operation on the basis of a symbolic interval base and degenerate interval exponent. Other constellations are not supported.

        Args:
            l (Interval): The base of the power operation.
            r (Interval): The exponent of the power operation.
            b (Boundary): The boundary ('Gültigkeitsbereich') of the power operation.

        Returns:
            _spow: A list of BoundedExpressions containing the result of the symbolic power operation (interval), possibly limited to subsets of the initial boundary.
        """
        if not self.is_constant(r.low):
            raise ValueError("The ExactExtremaEvaluator does not support symbolic exponents")

        if not (r.low == r.up):
            raise ValueError("The ExactExtremaEvaluator does not support non degenerate interval exponents")

        # 1. Determine Extrema
        _extrema = []

        pf = Power(
                Variable("x"),
                r.low
            )
        
        if isinstance(r.low, int) or r.low.is_integer():
            # 1. x^a, a > 0
            if r.low > 0:
                _extrema = [0]
            # 2. x^a, a = 0 => x^a = 1
            elif r.low == 0:
                return Interval(1)
            # 3. x^a, a < 0 = 1 / x^(|a|)
            else:
                return self._idiv(
                    1,
                    self._spow(
                        l=l,
                        r=Interval(abs(r.low)),
                        b=b
                    )
                )
        else:
            _extrema = self.__get_extrema(pf)

        return self.symbolic_expression(
            expr=pf,
            iv=l,
            extrema=_extrema,
            b=b
        )

    ####
    # UTILITY FUNCTIONS
    ####
    def _solve(self, f, b:Boundary=None) -> list:
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

    def _boundary_eval(self, f, bv: BoundedValue):
        """Evaluate a function for a given bounded value.

        Args:
            f: Function to evaluate.
            bv: Value at which the function is to be evaluated. This is either done as a limit, if the value is an open bound, or exactly.

        Returns:
            _boundary_eval: The evaluation of the function f for the bounded value.
        """
        from sympy.series.limits import limit
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

    def __merge_sections(self,existing:list,expansions:list) -> list:
        """Given a list composed of tuples specifying a boundary and the included extrema, this function merges this list with a list of new sections. Both lists have to specify sections over the same boundary."""
        # Merge with current res
        res = []

        for (exp_b, exp_in) in expansions:
            for (r_b,r_in) in existing:
                combined = r_b.intersect(exp_b)

                if combined is None:
                    continue
                    
                exp_in.extend(r_in)

                res.append(
                    (
                        combined,
                        exp_in
                    )
                )

        return res

    def __get_diff(self,pf:Expression):
        """Function to determine the derivative as sympy expression"""
        from vodes.symbolic.mapper.interop import ExactPymbolicToSympyMapper
        from pymbolic.mapper.differentiator import DifferentiationMapper as DM

        return ExactPymbolicToSympyMapper()(
            DM(Variable("x"))(pf)
        )

    def __get_extrema(self,pf:Expression) -> list:
        """Function to determine the (real) zero points of a given pymbolic Expression"""

        # 1. Differentiate the Operator
        f_diff = self.__get_diff(pf)
        self._logger.debug(f'Using {f_diff} for determining the extrema of the power operator.')

        # 2. Determine extrema
        # TODO : Evaluate constraint, if needed
        # TODO : e.g. x^3 -> differentiation between global/local extrema required
        # TODO : Ensure sorting
        return self._solve(f = f_diff)


