from abc import ABC, abstractmethod
from typing import List

from vodes.symbolic.translations.translation import Translation

# Custom Expression Library
from vodes.symbolic.expressions.bounded import BoundedExpression, Domain
from vodes.symbolic.expressions.interval import Interval

# Expression Library
from pymbolic.primitives import Expression

class ToScalar(Translation):
    """Converts expression to a scalar equivalent"""
    def translate(self, expr:BoundedExpression) -> List[BoundedExpression]:
        if isinstance(expr.expr,Interval):
            return [
                BoundedExpression(
                    expression = Interval(
                        lower=_translation(expr.expr.low,boundary=expr.bound,is_minimize=True),
                        upper=_translation(expr.expr.up,boundary=expr.bound,is_minimize=False)
                    ),
                    boundary=expr.bound
                )
            ]
        else:
            return [ 
                BoundedExpression(
                    expression=Interval(
                        lower=_translation(expr.expr,boundary=expr.bound,is_minimize=True),
                        upper=_translation(expr.expr,boundary=expr.bound,is_minimize=False)
                    ),
                    boundary=expr.bound
                )
            ]
        
def _translation(expr:Expression,boundary:Domain,is_minimize:bool) -> bool:
    from scipy.optimize import minimize_scalar
    from sympy import lambdify
    from vodes.symbolic.mapper.interop import ExactPymbolicToSympyMapper

    if not is_minimize:
        expr = expr * (-1)

    f = ExactPymbolicToSympyMapper()(expr)
    func = lambdify(f.free_symbols,f)
    
    # not rigorous, maybe return interval?
    # TODO : Evaluate
    res = minimize_scalar(func, bounds=(boundary.start, boundary.end), method='bounded')
    return func(res.x)

