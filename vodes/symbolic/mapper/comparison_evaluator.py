from typing import Dict, List, Tuple
from vodes.symbolic.analysis import Analysis, AnalysisConfig

# Assumption library
from vodes.symbolic.translations.to_scalar import ToScalar
from vodes.symbolic.properties.is_scalar import IsScalar
from vodes.symbolic.assumption import Assumption

# Custom Expression library
from vodes.symbolic.expressions.interval import Interval
from vodes.symbolic.expressions.bounded import BoundedVariable, BoundedExpression, Domain, DummyVariable, MachineError
from vodes.symbolic.expressions.trigonometric import sin,cos
from vodes.symbolic.expressions.infinity import Infinity, NegativeInfinity
from vodes.symbolic.utils import compare, le,ge,gt,lt,minimum,maximum

# Custom Mappers
from vodes.symbolic.mapper.interval_evaluator import IntervalEvaluator

# Expression Library
from pymbolic.primitives import Expression, Quotient, Variable, Power

def evaluate(expression, symbol:BoundedVariable, float:bool=False, context=None):
    if context is None:
        context = {}
    res = ComparisonEvaluator(context,symbol)(expression)

    # Two iterations of solver, if symbolic values are used for evaluation.
    # This allows to push the floating calculations further up.
    if float:
        return evaluate(res, symbol, float, context=context)
    else:
        return res

class ComparisonEvaluator(IntervalEvaluator):
    """Class for determining the exact boundaries of intervals on the basis of function analysis, powered by sympy."""

    def __init__(self, context: dict, symbol: BoundedVariable):
        super().__init__(context=context, symbol=symbol)
        self._assumptions["_ipow"] = (
                # Right Item (Exponent) has to be scalar value
                [],[
                    Assumption(
                        property=IsScalar()
                    )
                ]
            )

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
                request = []
                pfs = set()

                # We do this sorting, to allow extrema to be the first candidates in further methods.
                # This is advantageous for e.g. minima and maxima detection.
                # [-1,1] instead of [-1,-1](-1,1)[1,1] if boundaries are equal.
                # TODO : In the long run, this should be fixed within the minima and maxima detection (without adding further comparisons)

                # All extrema can be possible boundaries
                # e.g. [-x,x] ** 2 
                # != [(-x)**2,(x)**2]
                # == [min(0,(-x)**2,(x)**2),max(0,(-x)**2,(x)**2)] = [0,x**2]
                for extrema in contained:
                    extrema_expr = substitute(expr,{Variable("x"): extrema})

                    if extrema_expr in pfs:
                        continue

                    pfs.add(extrema_expr)
                    request.append(extrema_expr)

                for boundary in exprs:
                    if boundary in pfs:
                        continue

                    pfs.add(boundary)
                    request.append(boundary)


                self._logger.debug(f'Determining min/max for {list(map(str,request))}.')

                lower = self._minimum(request,subdomain)
                upper = self._maximum(request,subdomain)
                res.extend(
                    [ 
                        BoundedExpression(
                            expression=Interval(
                                lower=left,
                                upper=right
                            ),
                            boundary=boundary
                        ) for ((left,right),boundary) in merge(lower,upper)
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
        self._logger.info('==INTERVAL ADD==')
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
        self._logger.info('==INTERVAL SUB==')
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
        self._logger.info('==INTERVAL MUL==')

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
            ) for ((left,right),boundary) in merge(lower,upper)
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
        self._logger.info('==INTERVAL DIV==')
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
        self._logger.info('==INTERVAL POW==')

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
        from vodes.symbolic.utils import merge_unary
        self._logger.info('==INTERVAL ABS==')
        res = self._maximum(
            exprs=[
                i.up,
                (-1) * i.low
            ],
            boundary=d
        )

        # Maximum may return multiple equivalent expressions within an equivalent domain
        # Merge them using the domains
        return [
            BoundedExpression(
                expression=left,
                boundary=boundary
            ) for (left,boundary) in merge_unary(res)
        ] 

    def _inthroot(self, i:Interval, n:int, d:Domain) -> List[BoundedExpression]:
        from vodes.symbolic.expressions.infinity import NegativeInfinity
        self._logger.info('==INTERVAL NTHROOT==')

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

    def _isin(self, i:Interval, d:Domain) -> List[BoundedExpression]:
        self._logger.info('==INTERVAL SIN==')
        raise NotImplementedError()

    def _icos(self, i:Interval, d:Domain) -> List[BoundedExpression]:
        self._logger.info('==INTERVAL COS==')
        raise NotImplementedError()

    def _icontains(self, expr:Interval, val, d: Domain, incl:set=set(("up","low"))) -> list:
        """Determines if the provided symbolic interval contains the supplied value within the given domain."""
        from functools import cmp_to_key
        self._logger.info('==INTERVAL CONTAINS==')

        res = [
            # (
            #   Boundary ( start, end ),
            #   [ value ]
            # )
        ]

        if not ("up" in incl) :
            upper_contained = [d]
        else:
            # TODO : multivariate support
            upper = Analysis(expr.up,config=AnalysisConfig(d=d))
            upper_contained = upper.compare(val,lambda u,v: u>=v)

        if not ("low" in incl) :
            lower_contained = [d]
        else:
            # TODO : multivariate support
            lower = Analysis(expr.low,config=AnalysisConfig(d=d))
            lower_contained = lower.compare(val,lambda l,v: l<=v)

        for u in upper_contained:
            for l in lower_contained:
                combined = u.intersect(l)

                if combined is None:
                    continue
                
                res.append(
                    (combined,[val])
                )

        excluded = d.difference(list(map(lambda r:r[0],res)))
        for boundary in excluded:
            res.append(
                (boundary,[])
            )

        # Sort
        list.sort(res,key=cmp_to_key(lambda item1,item2: compare(item1[0],item2[0])))
                
        self._logger.debug(f'Determined sub-domains {res}')
        return res

    ####
    # SYMBOLIC EXPRESSION INTERFACE
    ####
    def _minimum(self, exprs:List[Expression],boundary:Domain) -> List[BoundedExpression]:
        exprs = list(set(exprs))

        (lassums) = self._assumptions.get("_minimum")
        candidates = [
            expr if not isinstance(expr, Interval) else expr.low for expr in self._apply_assumptions(lassums, exprs, boundary)
        ]

        self._logger.debug(f"Determining minima for {list(map(str,candidates))}") 
      
        # (1) Construct problem statements
        problems = [
            Analysis(pf=candidate,config=AnalysisConfig(d=boundary)) for candidate in candidates
            ]

        # (2) Analyze expression
        minima = self.analysis(problems=problems,d=boundary,method="min")
        
        # (3) Construct BoundedExpressions
        res = []
        
        for i in range(len(minima)):
            for boundary in minima[i]:
                res.append(
                    BoundedExpression(
                        expression=problems[i].expr,
                        boundary=boundary
                    )
                )

        self._logger.debug(f"MIN({list(map(str,exprs))}) => {list(map(str,res))}")
        return res

    def _maximum(self, exprs:List[Expression],boundary:Domain) -> List[BoundedExpression]:
        exprs = list(set(exprs))

        (lassums) = self._assumptions.get("_maximum")
        candidates = [
            expr if not isinstance(expr, Interval) else expr.low for expr in self._apply_assumptions(lassums, exprs, boundary)
        ]

        self._logger.debug(f"Determining maxima for {list(map(str,candidates))}") 
      
        # (1) Construct problem statements
        problems = [
            Analysis(pf=candidate,config=AnalysisConfig(d=boundary)) for candidate in candidates
            ]

        # (2) Analyze expression
        maxima = self.analysis(problems=problems,d=boundary,method="max")

        # (3) Construct BoundedExpressions
        res = []
        
        for i in range(len(maxima)):
            for boundary in maxima[i]:
                res.append(
                    BoundedExpression(
                        expression=problems[i].expr,
                        boundary=boundary
                    )
                )
        
        self._logger.debug(f"MAX({list(map(str,exprs))}) => {list(map(str,res))}")
        return res

    ####
    # INTERSECTION INTERFACE
    ####
    def analysis(self, problems: List[Expression], d:Domain, method:str="min"):
        """TODO"""
        from functools import reduce
        from vodes.symbolic.expressions.bounded import intersect

        if len(problems) == 0:
            raise ValueError("Supplied no candidates for analysis.")

        if method=="min":
            cmp = lambda l,r: l<=r
        elif method=="max":
            cmp = lambda l,r: l>=r
        else:
            raise ValueError(f"Not supporting method {method}")

        # Example : len(pr) = 4
        # 
        # | d                    pr[0] > pr[1]        pr[0] > pr[2]       pr[0] > pr[3] |
        # | ask(pr[0] <= pr[1])  d                    pr[1] > pr[2]       pr[1] > pr[3] |
        # | ask(pr[0] <= pr[2])  ask(pr[1] <= pr[2])  d                   pr[2] > pr[3] |
        # | ask(pr[0] <= pr[3])  ask(pr[1] <= pr[3])  ask(pr[2] <= pr[3]) d             |
        #

        res = [[[] for j in range(len(problems))] for i in range(len(problems))]

        extrema = [None for i in range(len(problems))]

        for i in range(len(problems)):
            # (1) Insert irrelevant comparisons 
            res[i][i] = [d]

            # (2) Fill upper triangular matrix incorporating previous results
            for j in range(i):
                res[i][j] = d.difference(res[j][i])

            # (3) Fill lower triangular matrix with comparison results
            for j in range(i+1,len(problems)):
                res[i][j] = problems[i].compare(problems[j],cmp,d=d)

            # (4) Store result
            extrema[i] = reduce(
                lambda a,b: intersect(a,b),
                res[i]
            )

        return extrema

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
        from vodes.symbolic.utils import eq

        rs_merged = [rs[0]]
        for i in range(1, len(rs)):
            if rs_merged[-1][0] == rs[i][0]:
                union = rs_merged[-1][1].union(rs[i][1])
                rs_merged[-1] = (
                    rs_merged[-1][0],
                    union[0]
                )

                for j in range(1,len(union)):
                    rs_merged.append(
                        (rs[i][0],union[j])
                    )
            else:  
                rs_merged.append(rs[i])

        return rs_merged

    def __merge_bounded_expressions(self, bexprs):
        """Merges bounds based on equivalency of expression"""
        return list(
            map(
                lambda tple: BoundedExpression(expression=tple[0],boundary=tple[1]),
                self.__merge_bounded_results(list(map(lambda bexpr : (bexpr.expr,bexpr.bound),bexprs)))
            ) 
        )