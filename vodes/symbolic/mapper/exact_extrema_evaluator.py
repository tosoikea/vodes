
from abc import abstractmethod, ABC

from vodes.symbolic.interval import Interval, ExactIntervalEvaluator
from vodes.symbolic.symbols import Boundary, BoundedValue, BoundedExpression
from vodes.symbolic.power import Power

# Expression library
from pymbolic.primitives import Variable, Expression, Quotient

class ExactExtremaEvaluator(ExactIntervalEvaluator, ABC):
    """Class for determining the exact boundaries of intervals on the basis of function analysis."""

    def __init__(self, context: dict, symbol: Variable):
        super().__init__(context=context, symbol=symbol)

            
    def general_symbolic_expression(self, expr:Expression, iv:Interval, b:Boundary) -> list:
        """
        TODO : Documentation ~ common_symbolic_expression. However, domains are obtained
        """
        from vodes.symbolic.mapper.interop import ExactPymbolicToSympyMapper as EPTS
        from sympy import S, symbols
        from sympy.calculus.util import continuous_domain
    
        sub_domains = []
        res = []

        # 0. Determine extrema
        extrema = self.__get_extrema(expr,b)

        # 1. Determine function bounds
        f = EPTS()(expr)
        d = continuous_domain(f,symbols("x"),S.Reals)
        self._logger.debug(f'Determined domain {d} for {f}')

        # 2. Determine undefined ranges
        invalid_bounds = self.__get_bounds(d.complement(S.Reals))

        # 3. Determine inclusions
        if len(invalid_bounds) == 0:
            sub_domains = [b]
        else:
            for invalid in invalid_bounds:
                self._logger.debug(f'Analyzing (partial) inclusion of {invalid} within {iv}')

                for (section_b, sections_vals) in self._icontains(iv=iv, xs=invalid, b=b):
                    if len(sections_vals) == 0:
                        sub_domains.append(section_b)
                    else:
                        self._logger.info(f'Dropping {section_b} as {f} contains values from {sections_vals} using {iv}')


        # 4. Evaluate expression
        for sub_domain in sub_domains:
            res.extend(
                self.common_symbolic_expression(
                    expr=expr,
                    iv=iv,
                    extrema=extrema,
                    b=sub_domain
                )
            )

        return res
            

    def common_symbolic_expression(self, expr:Expression, iv:Interval, extrema:list, b:Boundary) -> list:
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

        return self.__merge_bounded_expressions(res)

    def __get_sectioning(self,iv:Interval,extrema:list,b:Boundary) -> list:
        """Function to determine the sectioning based on when the interval contains the listed extremas."""
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

        for x in extrema:
            self._logger.debug(f'Analyzing zero point {x}.')
            res = self.__merge_sections(
                existing=res,
                expansions=self._icontains(
                    iv,
                    xs=x,
                    b=b
                )
            )

        #TODO : Minimize
        self._logger.debug(f'Determined the following extrema inclusions. {res}')
        return res

    

    ####
    # MIN/MAX INTERFACE
    ####
    @abstractmethod
    def _get_symbolic_interval(self, pfs:list, b: Boundary) -> list:
        """Construct symbolic intervals based on the included expressions within the given boundary. Possibly returns a list of multiple intervals split into sub-boundaries."""
        pass

    ####
    # SYMBOLIC INTERVAL INTERFACE
    ####
    def _sdiv(self, l:Interval, r:Interval, b: Boundary):
        """Interval division in the case of symbolic interval boundaries. In this case, the normal interval divison can not be determined, as knowledge about the inclusion of 0 is required.

        Args:
            l (Interval): The dividend of the division.
            r (Interval): The divisor of the division.
            b (Boundary): The boundary ('Gültigkeitsbereich') of the multiplication.

        Returns:
            _sdiv: A list of BoundedExpressions containing the result of the symbolic division (interval), possibly limited to subsets of the initial boundary.
        """
        res = []
        sections = self.__get_sectioning(iv=r, extrema=[0],b=b)

        for (b,xs) in sections:
            # Division undefined
            if 0 in xs:
                continue
            
            res.extend(
                self._imul(l, Interval(Quotient(1,r.up),Quotient(1,r.low)), b)
            )

        return res

    def _spow(self, l:Interval, r:Interval, b:Boundary):
        """Power operation on the basis of a symbolic interval base and degenerate interval exponent. Other constellations are not supported.

        Args:
            l (Interval): The base of the power operation.
            r (Interval): The exponent of the power operation.
            b (Boundary): The boundary ('Gültigkeitsbereich') of the power operation.

        Returns:
            _spow: A list of BoundedExpressions containing the result of the symbolic power operation (interval), possibly limited to subsets of the initial boundary.
        """
        from vodes.symbolic.utils import lt,gt

        if not self.is_constant(r.low):
            raise ValueError("The ExactExtremaEvaluator does not support symbolic exponents")

        if not (r.low == r.up):
            raise ValueError("The ExactExtremaEvaluator does not support non degenerate interval exponents")

        # Ensure scalar exponent
        exponent = r.low

        # x^a, a < 0 = 1 / x^(|a|)
        if lt(exponent,0):
            res = []
            bexprs = self._ipow(
                    l=l,
                    r=Interval(abs(exponent)),
                    b=b
                )
            
            for bexpr in bexprs:
                res.extend(
                    self._idiv(
                        Interval(1),
                        bexpr.expr,
                        bexpr.bound
                    )
                )
            
            return res

        pf = Power(
                Variable("x"),
                exponent
            )
            
        if isinstance(exponent, int):
            # 1. x^a, a > 0
            if gt(exponent,0):
                return self.common_symbolic_expression(
                    expr=pf,
                    iv=l,
                    extrema=[0],
                    # D(x^a) = R
                    b=b
                )
            # 2. x^a, a = 0 => x^a = 1
            else:
                return Interval(1)
        else:
            return self.general_symbolic_expression(
                expr=pf,
                iv=l,
                b=b
            )

    def _sinclusion(self, expr, val, incl, bf, b: Boundary) -> list:
        from vodes.symbolic.mapper.interop import ExactPymbolicToSympyMapper

        # 1. Determine sectioning
        # 
        # Every intersections between the expression and the desired value may indicate a change of inclusion.
        # We therefore start by splitting the boundary using these intersections and analyze the inclusion within them.
        # The inclusion is determined by evaluating the expression in the middle of the section and comparing it to the value.
        #
        sections = [
            #(start,end,included)
        ]

        f = ExactPymbolicToSympyMapper()(expr)
        
        # Ignore intersections etc. on boundary
        intersections = self._solve(
            f=f-val,
            b=Boundary(
                lower=BoundedValue(value=b.lower.value, open=True),
                upper=BoundedValue(value=b.upper.value, open=True)
            )
        )

        lower = b.lower.value

        # 1. Determine sectioning
        for x in intersections:
            if not b.contains(x):
                continue

            y = self._boundary_eval(
                f, 
                BoundedValue(
                    value=lower + (x-lower)/2,
                    open=False
                )
            )

            sections.append(
                (lower,x,incl(y, val))
            )
            lower = x

        if len(sections) == 0 or sections[-1][1] < b.upper.value:
            y = self._boundary_eval(
                f, 
                BoundedValue(
                    value=lower + (b.upper.value - lower)/2,
                    open=False
                )
            )

            sections.append(
                (lower, b.upper.value,incl(y, val))
            )

        
        self._logger.debug(f'Split {b} into {sections} for {f}')

        # 2. Determine sub-domains
        #
        # If we encounter two sections including the value, we collapse them.
        # The same does not apply to two sections not including the value. 
        # This is because they are seperated by an intersection between the expression and the value and therefore a sub-domain [val,val].

        res = [
            # (
            #   Boundary ( start, end ),
            #   [ value ]
            # )
        ]

        for i in range(len(sections)):
            (start,end,is_included) = sections[i]


            xs = [val] if is_included else []

            lower = BoundedValue(start, open=bf(is_included)) if i > 0 else b.lower
            upper = BoundedValue(end, open=bf(is_included)) if i < len(sections) - 1 else b.upper

            if i > 0:
                (_,_, prev_is_included) = sections[i-1]

                # Collapse
                if is_included and prev_is_included:
                    res[-1][0].upper = upper
                    continue

                # If previous was also not inclusive, append sub-domain with only intersection point as inclusive
                elif not is_included and not prev_is_included:
                    res.append(
                        (
                            Boundary(BoundedValue(start,False),BoundedValue(start,False)),
                            [val]
                        )
                    )
            
            res.append(
                (
                    Boundary(lower, upper),
                    xs
                )
            )
            
        self._logger.debug(f'Determined sub-domains {res}')

        return res

    ####
    # UTILITY FUNCTIONS
    ####
    def _purge_solutions(self, xs, b:Boundary=None) -> list:
        from sympy import im
        res = []

        for x in xs:
            #TODO : Make real/complex support a parameter choice
            if im(x) != 0:
                self._logger.debug(f'Dropping {x}.')
                continue

            if (b is None) or b.contains(x):
                res.append(x)

        return res

    def _map_solution(self, xs, b:Boundary=None) -> list:
        from sympy.sets.sets import Union
        from sympy.sets.fancysets import Reals

        res = [] 
        
        if xs.is_empty:
            self._logger.debug(f'No solution.')
            return res

        if xs.is_FiniteSet:
            return self._purge_solutions(xs, b=b)

        if isinstance(xs, Union):
            for set in xs.args:
                res.extend(
                    self._map_solution(set,b)
                )

            return res

        if isinstance(xs, Reals):
            self._logger.debug(f'No unique solution. Treating as if no solution is present.')
            return res

        self._logger.info(xs.__class__)
        raise ValueError("Found no solution.")

    def _solve(self, f, b:Boundary=None) -> list:
        from sympy.solvers import solveset
        from sympy import S
        self._logger.debug(f'Solving {f}')

        xs = solveset(f,domain=S.Reals)
        return self._map_solution(xs,b)

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

    def __get_bounded_value(self, value):
        from vodes.symbolic.symbols import Infinity, NegativeInfinity
        from sympy.core.numbers import Infinity as INF, NegativeInfinity as NINF
        from pymbolic.interop.sympy import SympyToPymbolicMapper
        

        if isinstance(value, INF):
            return Infinity()
        elif isinstance(value, NINF):
            return NegativeInfinity()
        else:
            return SympyToPymbolicMapper()(value)


    def __get_bounds(self, set) -> list:
        """Converts a sympy set to a boundary"""
        from sympy.sets.sets import EmptySet, Interval, Union

        res = []
        if isinstance(set,EmptySet):
            return res
        elif isinstance(set,Interval):
            res.append(
                Boundary(
                    lower=BoundedValue(self.__get_bounded_value(set.start),set.left_open),
                    upper=BoundedValue(self.__get_bounded_value(set.end),set.right_open)
                )
            )
        elif isinstance(set,Union):
            for x in set:
                res.extend(self.__get_bounds(x))

        return res


    def __get_extrema(self,pf:Expression,b:Boundary) -> list:
        """Function to determine the (real) zero points of a given pymbolic Expression"""
        from vodes.symbolic.mapper.interop import ExactPymbolicToSympyMapper as EPTS
        from pymbolic.primitives import Variable
        from pymbolic.mapper.differentiator import DifferentiationMapper as DM

        # 1. Differentiate the Operator
        pf_diff = DM(Variable("x"))(pf)
        f_diff = EPTS()(pf_diff)
        self._logger.debug(f'Using {f_diff} for determining the extrema of the power operator.')

        # 2. Determine extrema
        # TODO : Evaluate constraint, if needed
        # TODO : e.g. x^3 -> differentiation between global/local extrema required
        # TODO : Ensure sorting
        extrema = self._solve(f = f_diff)

        return extrema

