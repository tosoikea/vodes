from vodes.symbolic.properties.property import Property

# Custom Expression Library
from vodes.symbolic.expressions.interval import Interval
from vodes.symbolic.expressions.bounded import BoundedExpression

# Expression Library
from pymbolic.primitives import Expression

class IsScalar(Property):
    def verify(self, expr:BoundedExpression) -> bool:
        """Determines, if the expression is constant"""
        verify = [expr.expr.low,expr.expr.up] if isinstance(expr,Interval) else [expr.expr]

        return all(
            map(
                _verification,
                verify
            )
        )

def _verification(expr:Expression) -> bool:
    from pymbolic.primitives import is_constant
    from pymbolic.mapper.dependency import DependencyMapper

    return is_constant(expr) or not bool(DependencyMapper()(expr))