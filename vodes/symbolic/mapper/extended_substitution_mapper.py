from pymbolic.mapper.substitutor import SubstitutionMapper, make_subst_func

def substitute(expression, variable_assignments=None):
    if variable_assignments is None:
        variable_assignments = {}

    return ExtendedSubstitutionMapper(make_subst_func(variable_assignments))(expression)

class ExtendedSubstitutionMapper(SubstitutionMapper):
    ## FUNCTIONS
    def map_nthroot(self, expr):
        return expr.__class__(self.rec(expr.expr),expr.n)

    def map_sin(self, expr):
        return expr.__class__(self.rec(expr.expr))

    def map_cos(self, expr):
        return expr.__class__(self.rec(expr.expr))