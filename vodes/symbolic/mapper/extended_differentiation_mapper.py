from pymbolic.mapper.differentiator import DifferentiationMapper, map_math_functions_by_name
from pymbolic.primitives import Variable, Subscript, make_variable

def differentiate(expression,
                  variable,
                  allowed_nonsmoothness="none"):
    if not isinstance(variable, (Variable, Subscript)):
        variable = make_variable(variable)

    return ExtendedDifferentiationMapper(
        variable, map_math_functions_by_name, allowed_nonsmoothness=allowed_nonsmoothness
        )(expression)

class ExtendedDifferentiationMapper(DifferentiationMapper):
    """Class extending the default symbolic differentiation mapper to support custom expressions."""
    def map_sub(self,expr):
        from vodes.symbolic.expressions.primitives import Subtraction

        children = tuple([ self.rec(child) for child in expr.children ])

        return Subtraction(children)

    ## FUNCTIONS
    def __chain(self,outer,inner):
        from sympy import diff
        from vodes.symbolic.mapper.interop import ExactPymbolicToSympyMapper, ExactSympyToPymbolicMapper
        from vodes.symbolic.mapper.extended_substitution_mapper import substitute

        diff_outer = ExactSympyToPymbolicMapper()(
            diff(
                ExactPymbolicToSympyMapper()(outer)
            )
        )

        return substitute(
            diff_outer,
            {'x' : inner}
        ) * self.rec(inner)

    def map_nthroot(self, expr):
        from vodes.symbolic.expressions.nthroot import NthRoot

        return self.__chain(
            outer = NthRoot(Variable('x'),n=expr.n),
            inner = expr.expr
        )

    def map_sin(self, expr):
        from vodes.symbolic.expressions.trigonometric import cos
        
        return cos(expr.expr) * self.rec(expr.expr)

    def map_cos(self, expr):
        from vodes.symbolic.expressions.trigonometric import sin
        
        return (-sin(expr.expr)) * self.rec(expr.expr)