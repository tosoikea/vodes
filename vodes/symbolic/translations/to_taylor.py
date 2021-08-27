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
        from vodes.symbolic.utils import merge

        if isinstance(expr.expr,Interval):
            lower = [bexpr.expr.low for bexpr in self._translation(expr.expr.low,boundary=expr.bound)]
            upper = [bexpr.expr.up for bexpr in self._translation(expr.expr.up,boundary=expr.bound)]

            return [
                BoundedExpression(
                    expression=Interval(
                        lower=left,
                        upper=right
                    ),
                    boundary=boundary
                ) for (left,right,boundary) in merge(lower,upper)
            ]
        else:
            return self._translation(expr.expr,boundary=expr.bound)
        
    def _translation(self,expr:Expression,boundary:Domain) -> List[BoundedExpression]:
        from vodes.symbolic.analysis import Analysis, AnalysisConfig

        # TODO : Evaluate expansion point
        return Analysis(expr,config=AnalysisConfig(d=boundary)).taylor(n=self.n-1)

