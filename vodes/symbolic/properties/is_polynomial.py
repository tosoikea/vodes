from vodes.symbolic.properties.property import Property

# Custom Expression Library
from vodes.symbolic.expressions.interval import Interval
from vodes.symbolic.expressions.bounded import BoundedExpression

# Expression Library
from pymbolic.primitives import Expression

class IsPolynomial(Property):
    def __init__(self,n:int=4):
        assert n > 0
        self.n = n

    def verify(self, expr:BoundedExpression) -> bool:
        """Determines, if the expression is constant"""
        verify = [expr.expr.low,expr.expr.up] if isinstance(expr.expr,Interval) else [expr.expr]

        return all(
            map(
                self._verification,
                verify
            )
        )

    def _verification(self, expr:Expression) -> bool:
        from vodes.symbolic.analysis import Analysis
        return Analysis(expr).is_polynomial(n=self.n)