from typing import List

from vodes.symbolic.translations.translation import Translation

# Custom Expression Library
from vodes.symbolic.expressions.bounded import BoundedExpression, Domain
from vodes.symbolic.expressions.interval import Interval

# Expression Library
from pymbolic.primitives import Expression

class ToTaylor(Translation):
    def __init__(self,n:int=4,limit=2**16):
        assert n > 0
        self.n = n
        self.limit = limit
        
    """Provides taylor expansion for expression"""
    def translate(self, expr:BoundedExpression) -> List[BoundedExpression]:
        from vodes.symbolic.utils import merge

        if isinstance(expr.expr,Interval):
            lower = [bexpr.expr.low for bexpr in self._translation(expr.expr.low,boundary=expr.bound)]
            upper = [bexpr.expr.up for bexpr in self._translation(expr.expr.up,boundary=expr.bound)]
    
            res = [
                BoundedExpression(
                    expression=Interval(
                        lower=left,
                        upper=right
                    ),
                    boundary=boundary
                ) for ((left,right),boundary) in merge(lower,upper)
            ]
        else:
            res = self._translation(expr.expr,boundary=expr.bound)

        assert (len(res) == 1)

        return res
        
    def _translation(self,expr:Expression,boundary:Domain) -> List[BoundedExpression]:
        from vodes.symbolic.analysis import Analysis, AnalysisConfig

        # TODO : Evaluate expansion point
        res = Analysis(expr,config=AnalysisConfig(d=boundary,limit=self.limit)).taylor_model(n=self.n-1)
        return res

