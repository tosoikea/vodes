from typing import List, Tuple

# Assumption library
from vodes.symbolic.translations.to_scalar import ToScalar
from vodes.symbolic.properties.is_scalar import IsScalar
from vodes.symbolic.assumption import Assumption

# Custom Expression library
from vodes.symbolic.expressions.interval import Interval
from vodes.symbolic.expressions.bounded import BoundedVariable, BoundedExpression, Domain

# Custom Mappers
from vodes.symbolic.mapper.symbolic_interval_evaluator import ExactIntervalEvaluator

# Expression Library
from pymbolic.primitives import Expression, Quotient

class ScalarEvaluator(ExactIntervalEvaluator):
    """Class for determining the exact boundaries of intervals on the basis of function analysis."""

    def __init__(self, context: dict, symbol: BoundedVariable):
        super().__init__(context=context, symbol=symbol)
        self._assumptions = {
            "_minimum": [
                Assumption(
                    property=IsScalar(),
                    translation=ToScalar()
                )
            ],
            "_maximum": [
                Assumption(
                    property=IsScalar(),
                    translation=ToScalar()
                )
            ],
            "_idiv": [
                Assumption(
                    property=IsScalar(),
                    translation=ToScalar()
                )
            ],
            "_icontains":[
                Assumption(
                    property=IsScalar(),
                    translation=ToScalar()
                )
            ]
        }

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
        from vodes.symbolic.utils import le,ge
        rmin, rmax = r.low, r.up

        res = [
            BoundedExpression(expression=rmin,boundary=d),
            BoundedExpression(expression=rmax,boundary=d)
            ]
        for assumption in self._assumptions.get(self._idiv.__name__):
            res = assumption.validate(
                res
            )

        # 1. 0 \in r
        if le(rmin, 0) and ge(rmax, 0):
            # Infinity = Not Defined (for us)
            self._logger.warning(f'Dropping {r} for division, as it contains a 0.')
            return []
        # 2. 0 \not\in r
        else:
            return self._imul(l, Interval(Quotient(1,rmax),Quotient(1,rmin)), b)

    # TODO : Ensure that the function is defined on the domainTuple
    def _ipow(self, l:Interval, r:Interval, d: Domain) -> List[BoundedExpression]:
        raise NotImplementedError()
    
    def _iabs(self, i:Interval, d:Domain) -> List[BoundedExpression]:
        return self._maximum(
            exprs=[
                i.up,
                (-1) * i.low
            ],
            boundary=d
        )

    def _inthroot(self, i:Interval, n:int, d:Domain) -> List[BoundedExpression]:
        raise NotImplementedError()

    def _icontains(self, expr, val, incl, bf, d: Domain) -> list:
        """Determine the inclusion of any values of a boundary based on the inclusion of the upper and lower boundaries. For the interval to include a value, both boundaries have to include it."""  
        exprs = [
            BoundedExpression(expression=expr,boundary=d),
            ]
        for assumption in self._assumptions.get(self._icontains.__name__):
            exprs = assumption.validate(
                exprs
            )
        

        if len(exprs) != 1:
            raise ValueError(f"Unexpected {exprs} (len!=1)")

        if incl(expr.expr,val):
            return [
                (d,[val])
            ]
        else:
            return [
                (d,[])
            ]

    

    ####
    # SYMBOLIC EXPRESSION INTERFACE
    ####
    def _minimum(self, exprs:List[Expression],boundary:Domain) -> List[BoundedExpression]:
        from vodes.symbolic.utils import minimum

        res = [BoundedExpression(expression=expr,boundary=boundary) for expr in exprs]
        for assumption in self._assumptions.get(self._minimum.__name__):
            res = assumption.validate(
                res
            )

        vals = [
            val.expr if not isinstance(val.expr,Interval) else val.expr.low for val in res
        ]         
        res = minimum(vals)

        self._logger.debug(f'Minimum {res} determined from {exprs}')
        return [
            BoundedExpression(
                expression=res,
                boundary=boundary
            )
        ]

    def _maximum(self, exprs:List[Expression],boundary:Domain) -> List[BoundedExpression]:
        from vodes.symbolic.utils import maximum

        res = [BoundedExpression(expression=expr,boundary=boundary) for expr in exprs]
        for assumption in self._assumptions.get(self._maximum.__name__):
            res = assumption.validate(
                res
            )

        vals = [
            val.expr if not isinstance(val.expr,Interval) else val.expr.up for val in res
        ]
        res = maximum(vals)

        self._logger.debug(f'Maximum {res} determined from {exprs}')
        return [
            BoundedExpression(
                expression=res,
                boundary=boundary
            )
        ]

    