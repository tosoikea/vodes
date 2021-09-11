import logging
from typing import List, Tuple
from vodes.symbolic.expressions.bounded import BoundedExpression, BoundedVariable, Domain, DummyVariable
from vodes.symbolic.mapper.interval_evaluator import IntervalEvaluator
from vodes.symbolic.expressions.infinity import NegativeInfinity
from vodes.symbolic.translations.nop_translation import NOPTranslation

# Assumption library
from vodes.symbolic.translations.to_scalar import ToScalar
from vodes.symbolic.properties.is_scalar import IsScalar
from vodes.symbolic.assumption import Assumption

# Custom Expression library
from vodes.symbolic.expressions.interval import Interval
from vodes.symbolic.expressions.nthroot import NthRoot

# Custom Mappers

# Expression Library
from pymbolic.primitives import Expression, Quotient

# Expression Mapper

def evaluate(expression, context=None, float:bool=False):
    if context is None:
        context = {}
    res = ScalarEvaluator(context,DummyVariable())(expression)

    # Two iterations of solver, if symbolic values are used for evaluation.
    # This allows to push the floating calculations further up.
    if float:
        return evaluate(res, context=context)
    else:
        return res

class ScalarEvaluator(IntervalEvaluator):
    """Class for determining the exact boundaries of intervals on the basis of function analysis."""
    def __init__(self, context: dict, symbol:BoundedVariable):
        super().__init__(context=context, symbol=symbol)
        self._assumptions["_minimum"] = ([
                Assumption(
                    property=IsScalar(),
                    translation=ToScalar()
                )
            ]) 

        self._assumptions["_maximum"] = ([
            Assumption(
                property=IsScalar(),
                translation=ToScalar()
            )
        ])

        self._assumptions["_icontains"] = ([
            Assumption(
                property=IsScalar(),
                translation=ToScalar()
            )
        ])

        self._logger = logging.getLogger(__name__)
        self._context = context

    ####
    # INTERVAL INTERFACE
    ####
    def _iadd(self, l:Interval, r:Interval, d:Domain) -> List[BoundedExpression]:
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

    def _isub(self, l:Interval, r:Interval, d:Domain) -> List[BoundedExpression]:
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

    def _idiv(self, l:Interval, r:Interval, d:Domain) -> List[BoundedExpression]:
        """Interval division as defined within https://epubs.siam.org/doi/abs/10.1137/1.9780898717716.ch2.
        
        Args:
            l (Interval): The dividend of the division.
            r (Interval): The divisor of the division.
            b (Boundary): The boundary ('Gültigkeitsbereich') of the division.

        Returns:
            _idiv: A list of BoundedExpressions containing the result of the division (interval), possibly limited to subsets of the initial boundary, that may in sum not be equal to the initial boundary.
        """
        rmin, rmax = r.low, r.up

        # 1. 0 \in r
        assoc = self.contains(r, 0 , d)
        if len(assoc) > 1 or len(assoc[0][1]) > 0:
            # Infinity = Not Defined (for us)
            self._logger.warning(f'Dropping {r} for division, as it contains a 0.')
            return []
        # 2. 0 \not\in r
        else:
            return self._imul(l, Interval(Quotient(1,rmax),Quotient(1,rmin)),d)

    def _imul(self, l:Interval, r:Interval, d: Domain) -> List[BoundedExpression]:
        """Interval multiplication as defined within https://epubs.siam.org/doi/abs/10.1137/1.9780898717716.ch2. Because symbolic intervals are supported, the inheriting classes define custom evaluations for symbolic cases.
        
        Args:
            l (Interval): The left parameter of the multiplication.
            r (Interval): The right parameter of the multiplication.
            b (Boundary): The boundary ('Gültigkeitsbereich') of the multiplication.

        Returns:
            _imul: A list of BoundedExpressions containing the result of the multiplication (interval), possibly limited to subsets of the initial boundary.
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
            ) for ((left,right),boundary) in merge(lower,upper)
        ]

    def _ipow(self, l:Interval, r:Interval, d:Domain) -> List[BoundedExpression]:
        from vodes.symbolic.mapper.interop import ExactPymbolicToSympyMapper, ExactSympyToPymbolicMapper

        if not (r.low) == (r.up):
            raise ValueError("Not supporting non degenerate exponents")

        return [
            BoundedExpression(
                expression=ExactSympyToPymbolicMapper()(
                    ExactPymbolicToSympyMapper()(l) ** ExactPymbolicToSympyMapper()(r)
                ),
                boundary=d
            )
        ]
    
    def _iabs(self, i:Interval, d:Domain) -> List[BoundedExpression]:
        return self._maximum(
            exprs=[
                i.up,
                (-1) * i.low
            ],
            d=d
        )

    def _inthroot(self, i:Interval, n:int, d:Domain) -> Expression:
        # strictly monotonic
        # TODO : Assert positive value within i
        return [
            BoundedExpression(
                expression = Interval(
                    NthRoot(i.low,n),
                    NthRoot(i.up,n)
                ),
                boundary=d
            )
        ]
        
    def _isin(self, i:Interval, d:Domain) -> Expression:
        pass

    def _icos(self, i:Interval, d:Domain) -> Expression:
        pass 
    
    def _icontains(self, expr:Interval, val, d: Domain, incl:set=set(("up","low"))) -> list:
        """Determine the inclusion of any values of a boundary based on the inclusion of the upper and lower boundaries. For the interval to include a value, both boundaries have to include it."""  
        from vodes.symbolic.utils import le,ge

        if ("up" in incl) and not ge(expr.up,val):
            return [(d,[])]
        elif ("low" in incl) and not le(expr.low,val):
            return [(d,[])]
        else:
            return [(d,[val])]

    ####
    # SYMBOLIC EXPRESSION INTERFACE
    ####
    def _minimum(self, exprs:List[Expression],d:Domain) -> List[BoundedExpression]:
        from vodes.symbolic.utils import minimum

        res = minimum(exprs)

        self._logger.debug(f'Minimum {res} determined from {list(map(str,exprs))}')
        return [
            BoundedExpression(
                expression=res,
                boundary=d
            )
        ]

    def _maximum(self, exprs:List[Expression],d:Domain) -> List[BoundedExpression]:
        from vodes.symbolic.utils import maximum

        res = maximum(exprs)

        self._logger.debug(f'Maximum {res} determined from {list(map(str,exprs))}')
        return [
            BoundedExpression(
                expression=res,
                boundary=d
            )
        ]

    