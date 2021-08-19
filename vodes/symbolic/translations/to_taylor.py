from abc import ABC, abstractmethod
from typing import List

from vodes.symbolic.translations.translation import Translation

# Custom Expression Library
from vodes.symbolic.expressions.bounded import BoundedExpression, Domain
from vodes.symbolic.expressions.interval import Interval

# Expression Library
from pymbolic.primitives import Expression

class ToTaylor(Translation):
    def __init__(self,n:int=4):
        assert n > 0
        self.n = n
        
    """Provides taylor expansion for expression"""
    def translate(self, expr:BoundedExpression) -> List[BoundedExpression]:
        if isinstance(expr.expr,Interval):
            return [
                BoundedExpression(
                    expression = Interval(
                        lower=self._translation(expr.expr.low,boundary=expr.bound,is_minimize=True),
                        upper=self._translation(expr.expr.up,boundary=expr.bound,is_minimize=False)
                    ),
                    boundary=expr.bound
                )
            ]
        else:
            return [ 
                BoundedExpression(
                    expression=Interval(
                        lower=self._translation(expr.expr,boundary=expr.bound,is_minimize=True),
                        upper=self._translation(expr.expr,boundary=expr.bound,is_minimize=False)
                    ),
                    boundary=expr.bound
                )
            ]
        
    def _translation(self,expr:Expression,boundary:Domain,is_minimize:bool) -> bool:
        from vodes.symbolic.analysis import Analysis
        from pymbolic.primitives import Quotient
        from fractions import Fraction
        from pymbolic import evaluate

        # Prevent Floating Number -> closes rational
        midf = Fraction(evaluate(boundary.start + (boundary.end - boundary.start)/2))#.limit_denominator(...)
        mid = Quotient(midf.numerator, midf.denominator)

        mid = 0
        taylor_expr = Analysis(expr).taylor(a=mid,n=2)
        return taylor_expr

