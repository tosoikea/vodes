from pymbolic.mapper.substitutor import SubstitutionMapper, make_subst_func

def substitute(expression, variable_assignments=None):
    if variable_assignments is None:
        variable_assignments = {}

    return ExtendedSubstitutionMapper(make_subst_func(variable_assignments))(expression)

class ExtendedSubstitutionMapper(SubstitutionMapper):
    def map_variable(self, expr):
        result = self.subst_func(expr)
        if result is not None:
            return result
        else:
            return expr

    def map_sub(self,expr):
        from vodes.symbolic.expressions.primitives import Subtraction

        children = tuple([ self.rec(child) for child in expr.children ])

        return Subtraction(children)

    def map_interval(self,expr):
        from vodes.symbolic.expressions.interval import Interval
        
        return Interval(
            self.rec(expr.low),
            self.rec(expr.up)
        )

    ## FUNCTIONS
    def map_nthroot(self, expr):
        return expr.__class__(self.rec(expr.expr),expr.n)

    def map_sin(self, expr):
        return expr.__class__(self.rec(expr.expr))

    def map_cos(self, expr):
        return expr.__class__(self.rec(expr.expr))