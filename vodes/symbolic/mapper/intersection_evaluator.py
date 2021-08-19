from typing import List, Tuple
from vodes.symbolic.analysis import Analysis

# Assumption library
from vodes.symbolic.translations.to_scalar import ToScalar
from vodes.symbolic.properties.is_scalar import IsScalar
from vodes.symbolic.assumption import Assumption

# Custom Expression library
from vodes.symbolic.expressions.interval import Interval
from vodes.symbolic.expressions.bounded import BoundedVariable, BoundedExpression, Domain
from vodes.symbolic.utils import le,ge,gt,lt,minimum,maximum

# Custom Mappers
from vodes.symbolic.mapper.symbolic_interval_evaluator import ExactIntervalEvaluator

# Expression Library
from pymbolic.primitives import Expression, Quotient, Variable, Power

class IntersectionEvaluator(ExactIntervalEvaluator):
    """Class for determining the exact boundaries of intervals on the basis of function analysis.
    TODO : Describe in depth, firstly in thesis."""

    def __init__(self, context: dict, symbol: BoundedVariable):
        super().__init__(context=context, symbol=symbol)
        self._assumptions = {
            "_minimum": [],
            "_maximum": [],
            "_ipow": [
                Assumption(
                    property=IsScalar()
                )
            ],
            "_icontains":[]
        }

    def common_symbolic_expression(self, expr:Expression, iv:Interval, d:Domain, extrema:list=[], invalid_bounds:List[Domain]=None) -> list:
        """
        TODO : Documentation
        """
        from pymbolic import substitute
        from vodes.symbolic.utils import merge

        res = []
        sub_domains = []

        if invalid_bounds is None or len(invalid_bounds) == 0:
            sub_domains = [d]
        else:
            for invalid in invalid_bounds:
                self._logger.debug(f'Analyzing (partial) inclusion of {invalid} within {iv}')

                for (section_b, sections_vals) in self.contains(iv=iv,xs=invalid,d=d):
                    if len(sections_vals) == 0:
                        sub_domains.append(section_b)
                    else:
                        self._logger.info(f'Dropping {section_b} as {expr} contains values from {sections_vals} using {iv}')

        for sub_domain in sub_domains:
            # 1. Determine Sections
            sections = self.__get_sectioning(iv=iv, extrema=extrema,d=sub_domain)

            # 2. Determine min, max for sections
            exprs = [
                substitute(expr,{Variable("x"): iv.low}),
                substitute(expr,{Variable("x"): iv.up})
            ]

            for (subdomain,contained) in sections:
                pfs = set(exprs)

                # All extrema can be possible boundaries
                # e.g. [-x,x] ** 2 
                # != [(-x)**2,(x)**2]
                # == [min(0,(-x)**2,(x)**2),max(0,(-x)**2,(x)**2)] = [0,x**2]
                for extrema in contained:
                    pfs.add(substitute(expr,{Variable("x"): extrema}))

                self._logger.debug(f'Determining min/max for {pfs}.')

                lower = self._minimum(list(pfs),subdomain)
                upper = self._maximum(list(pfs),subdomain)

                res.extend(
                    [ 
                        BoundedExpression(
                            expression=Interval(
                                lower=left,
                                upper=right
                            ),
                            boundary=boundary
                        ) for (left,right,boundary) in merge(lower,upper)
                    ]
                )

        return self.__merge_bounded_expressions(res) 

    ####
    # INTERVAL INTERFACE
    ####
    def _iadd(self, l:Interval, r:Interval, d: Domain) -> List[BoundedExpression]:
        """Interval addition as defined within https://epubs.siam.org/doi/abs/10.1137/1.9780898717716.ch2.
        
        Args:
            l (Interval): The left parameter of the addition.
            r (Interval): The right parameter of the addition.
            b (Boundary): The boundary ('Gültigkeitsbereich') of the addition.

        Returns:
            _iadd: A list of BoundedExpressions containing the result of the addition (interval), possibly limited to subsets of the initial boundary, that may in sum not be equal to the initial boundary.
        """
        return [
            BoundedExpression(
                expression=Interval(l.low + r.low, l.up + r.up),
                boundary=d
            )
        ]

    def _isub(self, l:Interval, r:Interval, d: Domain) -> List[BoundedExpression]:
        """Interval substitution as defined within https://epubs.siam.org/doi/abs/10.1137/1.9780898717716.ch2.
        
        Args:
            l (Interval): The left parameter of the substitution.
            r (Interval): The right parameter of the substitution.
            b (Boundary): The boundary ('Gültigkeitsbereich') of the substitution.

        Returns:
            _isub: A list of BoundedExpressions containing the result of the substitution (interval), possibly limited to subsets of the initial boundary, that may in sum not be equal to the initial boundary.
        """
        return [
            BoundedExpression(
                expression=Interval(l.low - r.up, l.up - r.low),
                boundary=d
            )
        ]

    def _imul(self, l:Interval, r:Interval, d: Domain) -> List[BoundedExpression]:
        """Interval multiplication as defined within https://epubs.siam.org/doi/abs/10.1137/1.9780898717716.ch2. Because symbolic intervals are supported, the inheriting classes define custom evaluations for symbolic cases.
        
        Args:
            l (Interval): The left parameter of the multiplication.
            r (Interval): The right parameter of the multiplication.
            b (Boundary): The boundary ('Gültigkeitsbereich') of the multiplication.

        Returns:
            _imul: A list of BoundedExpressions containing the result of the symbolic multiplication (interval), possibly limited to subsets of the initial boundary.
        """
        from vodes.symbolic.utils import merge

        exprs = [
                l.low * r.low, 
                l.low * r.up, 
                l.up * r.low, 
                l.up * r.up
            ]

        lower = self._minimum(exprs,d)
        upper = self._maximum(exprs,d)

        return [
            BoundedExpression(
                expression=Interval(
                    lower=left,
                    upper=right
                ),
                boundary=boundary
            ) for (left,right,boundary) in merge(lower,upper)
        ]

    def _idiv(self, l:Interval, r:Interval, d: Domain) -> List[BoundedExpression]:
        """Interval division as defined within https://epubs.siam.org/doi/abs/10.1137/1.9780898717716.ch2.
        
        Args:
            l (Interval): The dividend of the division.
            r (Interval): The divisor of the division.
            b (Boundary): The boundary ('Gültigkeitsbereich') of the division.

        Returns:
            _idiv: A list of BoundedExpressions containing the result of the division (interval), possibly limited to subsets of the initial boundary, that may in sum not be equal to the initial boundary.
        """
        rmin, rmax = r.low, r.up
        inclusion = self.contains(iv=r,xs=0,d=d)

        res = []

        for (sub_domain,included) in inclusion:
            if 0 in included:
                self._logger.warning(f'Dropping sub-domain {sub_domain} for division with {r}, as it contains a 0.')
                continue

            res.extend(
                self._imul(l, Interval(Quotient(1,rmax),Quotient(1,rmin)), sub_domain)
            )

        return res

    def _ipow(self, l:Interval, r:Interval, d: Domain) -> List[BoundedExpression]:
        """Power operation on the basis of a symbolic interval base and degenerate interval exponent. Other constellations are not supported.

        Args:
            l (Interval): The base of the power operation.
            r (Interval): The exponent of the power operation.
            b (Boundary): The boundary ('Gültigkeitsbereich') of the power operation.

        Returns:
            _ipow: A list of BoundedExpressions containing the result of the symbolic power operation (interval), possibly limited to subsets of the initial boundary.
        """
        from vodes.symbolic.utils import lt,gt,eq
        from vodes.symbolic.mapper.simplification_mapper import simplify

        # TODO : Maybe convert to scalar interval, requires support for interval exponentiation
        exprs = [
            BoundedExpression(expression=r,boundary=d),
        ]
        for assumption in self._assumptions.get(self._icontains.__name__):
            exprs = assumption.validate(
                exprs
            )

        if not (r.low == r.up):
            raise ValueError("The IntersectionEvaluator does not support non degenerate interval exponents")

        res = []
        exponent = simplify(r.low)

        if lt(exponent,0):
            bexprs = self._ipow(
                    l=l,
                    r=Interval(abs(exponent)),
                    d=d
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

        # TODO : Simplify
        # Reason for limiting b^x to x \in N :
        # 
        # The definitions for x \in R introduce further constraints, that are best considered by using the underlying constructs.
        # e.g. b^r == exp(r * ln(b)), b \in R^+
        pf = Power(
            Variable("x"),
            exponent
        )

        # 1. b^x, x \in N
        if isinstance(exponent, int):
            # a.) b^x, x > 0
            if gt(exponent,0):
                extrema = [0] if eq(exponent % 2,0) else []
                return self.common_symbolic_expression(
                    expr=pf,
                    iv=l,
                    extrema=extrema,
                    d=d
                )
            # b.) b^x, x = 0 => b^x = 1
            else:
                return Interval(1)
        else:
            raise ValueError(f"Not supporting exponent {exponent}({exponent.__class__}). Please use the underyling constructs (exp,log,nth-root)...")

    def _iabs(self, i:Interval, d:Domain) -> List[BoundedExpression]:
        return self._maximum(
            exprs=[
                i.up,
                (-1) * i.low
            ],
            boundary=d
        )

    def _inthroot(self, i:Interval, n:int, d:Domain) -> List[BoundedExpression]:
        from vodes.symbolic.expressions.infinity import NegativeInfinity

        invalid_bounds = Domain(
                        start=NegativeInfinity(),
                        end=0,
                        left_open=True,
                        right_open=True
                    ) if n % 2 == 0 else None

        inclusion = self.contains(iv=i,xs=invalid_bounds,d=d)

        res = []

        for (sub_domain,included) in inclusion:
            if len(included) > 0:
                self._logger.warning(f'Dropping sub-domain {sub_domain} for {n}th-root with {i}, as it is undefined.')
                continue
            
            # nth-root is monotonic increasing
            # TODO : Verify this assumption
            res.append(
                BoundedExpression(
                    expression=Interval(
                        Power(base=i.low,exponent=Quotient(1,n)),
                        Power(base=i.up,exponent=Quotient(1,n))
                    ),
                    boundary=sub_domain
                )
            )

        return res
        

    def _icontains(self, expr:Expression, val, incl, bf, d: Domain) -> list:
        """Determines if the provided symbolic boundary expr contains the supplied value within the given domain using the supplied inclusion function."""
        exprs = [
            BoundedExpression(expression=expr,boundary=d),
            ]
        for assumption in self._assumptions.get(self._icontains.__name__):
            exprs = assumption.validate(
                exprs
            )

        res = [
            # (
            #   Boundary ( start, end ),
            #   [ value ]
            # )
        ]

        for bexpr in exprs:
            # 1. Determine sectioning
            # 
            # Every intersections between the expression and the desired value may indicate a change of inclusion.
            # We therefore start by splitting the boundary using these intersections and analyze the inclusion within them.
            # The inclusion is determined by evaluating the expression in the middle of the section and comparing it to the value.
            #
            sections = [
                #(start,end,included)
            ]

            problem = Analysis(bexpr.expr)
            
            domains = problem.equals(
                val,
                # Ignore intersections etc. on boundary
                Domain(
                    start=d.start,
                    end=d.end,
                    left_open=True,
                    right_open=True
                )
            )

            lower = d.start

            # 1. Determine sectioning
            for domain in domains:
                y = problem.evaluate(
                    lower + (domain.start-lower)/2
                )

                # TODO : Evaluate domain usage, especially for open doamins
                sections.append(
                    (lower,domain.start,incl(y, val))
                )
                sections.append(
                    (domain.start,domain.end,True)
                )

                lower = domain.end

            if len(sections) == 0 or le(sections[-1][1],d.end):
                y = problem.evaluate(
                    lower + (d.end-lower)/2
                )

                sections.append(
                    (lower,d.end,incl(y, val))
                )

            
            self._logger.debug(f'Split {d} into {sections} for {expr}')

            # 2. Determine sub-domains
            #
            # If we encounter two sections including the value, we collapse them.
            # The same does not apply to two sections not including the value. 
            # This is because they are seperated by an intersection between the expression and the value and therefore a sub-domain [val,val].
            for i in range(len(sections)):
                (start,end,is_included) = sections[i]


                xs = [val] if is_included else []

                (lv,lo) = (start, bf(is_included)) if i > 0 else (d.start, d.left_open)
                (rv,ro) = (end, bf(is_included)) if i < len(sections) - 1 else (d.end, d.right_open)

                if i > 0:
                    (_,_, prev_is_included) = sections[i-1]

                    # Collapse
                    if is_included and prev_is_included:
                        res[-1][0].end = rv
                        res[-1][0].right_open = ro
                        continue

                res.append(
                    (
                        Domain(
                            start=lv,
                            end=rv,
                            left_open=lo,
                            right_open=ro
                        ),
                        xs
                    )
                )
            
        self._logger.debug(f'Determined sub-domains {res}')
        return res

    ####
    # SYMBOLIC EXPRESSION INTERFACE
    ####
    def _minimum(self, exprs:List[Expression],boundary:Domain) -> List[BoundedExpression]:
        res = [BoundedExpression(expression=expr,boundary=boundary) for expr in exprs]
        for assumption in self._assumptions.get(self._minimum.__name__):
            res = assumption.validate(
                res
            )

        vals = [
            val.expr if not isinstance(val.expr,Interval) else val.expr.low for val in res
        ] 

        self._logger.debug(f"Determining minimum for {vals}") 
      
        # (1) Construct problem statements
        problems = [Analysis(pf=val) for val in vals]

        # (2) Determine intersection for possible min/max switch
        intersections = self.__intersections(problems,boundary)
        self._logger.debug(f'Intersections : {intersections}')

        # (3) Analyze problems using intersections
        res = [
           BoundedExpression(
               expression=problems[extremum].expr,
               boundary=subdomain
           )
           for (extremum, subdomain) in self.__analysis(self.__min_eval,intersections=intersections,problems=problems,d=boundary)
        ]

        self._logger.debug(f"MIN({exprs}) => {res}")
        return res

    def _maximum(self, exprs:List[Expression],boundary:Domain) -> List[BoundedExpression]:
        res = [BoundedExpression(expression=expr,boundary=boundary) for expr in exprs]
        for assumption in self._assumptions.get(self._maximum.__name__):
            res = assumption.validate(
                res
            )

        vals = [
            val.expr if not isinstance(val.expr,Interval) else val.expr.low for val in res
        ]  

        self._logger.debug(f"Determining maximum for {vals}")      
      
        # (1) Construct problem statements
        problems = [Analysis(pf=val) for val in vals]

        # (2) Determine intersection for possible min/max switch
        intersections = self.__intersections(problems,boundary)
        self._logger.debug(f'Intersections : {intersections}')

        # (3) Analyze problems using intersections
        res = [
           BoundedExpression(
               expression=problems[extremum].expr,
               boundary=subdomain
           )
           for (extremum, subdomain) in self.__analysis(self.__max_eval,intersections=intersections,problems=problems,d=boundary)
        ]

        self._logger.debug(f"MAX({list(map(str,exprs))}) => {list(map(str,res))}")
        return res

    ####
    # INTERSECTION INTERFACE
    ####
    def __min_eval(self, x0, x1,fs:List[Analysis]):
        """Evaluate the minimum function within the given boundary."""
        return self.__eval(fs=fs,x0=x0, x1=x1, order=lambda a, b: lt(a, b))

    def __max_eval(self, x0, x1,fs:List[Analysis]):
        """Evaluate the maximum function within the given boundary."""
        return self.__eval(fs=fs,x0=x0, x1=x1, order=lambda a, b: gt(a, b))

    def __eval(self, order, x0, x1, fs:List[Analysis]):
        """Evaluate problems between two values and determine the extremum using the given order function.

        Args:
            fs: Problems to evaluate.
            x0: Start of the section
            x1: End of the section
            order: Defines the comparison between function value (e.g. <, >).
        
        Returns:
            __eval: 
        """
        values = [None] * len(fs)
        x = x0 + (x1-x0)/2

        self._logger.debug(f'Evaluating {list(map(lambda f:f.func,fs))} at {x} ({x0}, {x1})')

        for i in range(0, len(fs)):
            values[i] = fs[i].evaluate(x)
        
        result = (0, values[0])
        for i in range(1, len(fs)):
            if order(values[i], result[1]):
                result = (i, values[i])
            elif result[1] == values[i]:
                # TODO : Take this from the intersection, where equivalence would be returned in terms of domain -> Less function evaluation
                # As the evaluation takes places between two intersections, the functions must be equivalent, as only these cases of identical values are not contained within the intersections.
                self._logger.info(f'Encountered equivalent functions {fs[result[0]].func} and {fs[i].func} within the section.')

        return result[0]

    def __analysis(self, eval, intersections: list, problems: List[Analysis], d:Domain):
        """Evaluate the extremas (e.g. minima, maxima) of the given functions within the boundary and return a list of them with their respective boundary."""
        res = []

        lv,lo = d.start,d.left_open
        candidates = [i for i in range(len(problems))]
        
        # There is still domain left to analyze
        is_analyzed = False
        while not is_analyzed:
            x0 = lv

            # Determine the next possible intersection
            xs = [item for sublist in intersections for subsublist in sublist for item in subsublist if gt(item,x0)]
            xs.append(d.end)

            x1 = minimum(xs)

            # Determine extremum
            extr_i = candidates[
                eval(
                    x0=x0,
                    x1=x1,
                    fs=[problems[i] for i in candidates]
                    )
                ]

            # Determine sub-domain
            rv,ro = d.end, d.right_open
            for i in range(len(intersections[extr_i])):
                xs = [item for item in intersections[extr_i][i] if gt(item,x0) and le(item,rv)]

                if len(xs) == 0:
                    continue
                else:
                    t = minimum(xs)

                    if lt(t, rv):
                        # mawarny be tangent
                        candidates = [extr_i, i]

                        rv,ro = t,True
                    else:
                        candidates.append(i)
                    
            boundary=Domain(lv,rv,lo,ro) 
            self._logger.debug(f'Extrema : {problems[extr_i].func} within {boundary}')

            res.append(
                (extr_i, boundary)
            )

            lv,lo = rv, not ro

            # Finalize, if domain was fully analyzed
            is_analyzed = lv == d.end

        return res

    def __intersection(self, f:Analysis, g:Analysis, d: Domain):
        """Determine the real intersections of two functions for the given boundary. Intersections at the boundary values are ignored.

        Args:
            f: Function to intersect with g.
            g: Function to intersect with f.
            b: Boundary to limit the possible values for the intersections.

        Returns:
            __intersection: A (possibly empty) list of intersections between f and g.
        """
        answer = f.equals(
            g.func,
            Domain(
                d.start,
                d.end,
                True,
                True
            )
        )

        res = []

        # TODO : Make use of domains
        for domain in answer:
            if domain.start != domain.end:
                self._logger.info(f'Dropping {domain} as domains are not yet fully supported.')
                continue

            res.append(domain.start)
        
        return res

    def __intersections(self, fs:List[Analysis], d: Domain):
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

                intersections[i][j] = self.__intersection(fs[i], fs[j], d)
                intersections[j][i] = intersections[i][j]
        
        return intersections


    ####
    # UTILITY FUNCTIONS
    ####
    def __get_sectioning(self,iv:Interval,extrema:list,d:Domain) -> list:
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
                d,
                []
            )
        ]

        for x in extrema:
            self._logger.debug(f'Analyzing zero point {x}.')
            res = self.__merge_sections(
                existing=res,
                expansions=self.contains(
                    iv,
                    xs=x,
                    d=d
                )
            )

        #TODO : Minimize
        self._logger.debug(f'Determined the following extrema inclusions. {res}')
        return res

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


    