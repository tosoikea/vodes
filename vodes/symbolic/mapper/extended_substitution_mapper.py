from pymbolic.mapper.substitutor import SubstitutionMapper

class ExtendedSubstitutionMapper(SubstitutionMapper):
    ## FUNCTIONS
    def map_nthroot(self, expr):
        return expr.__class__(self.rec(expr.expr),expr.n)

    def map_sin(self, expr):
        return expr.__class__(self.rec(expr.expr))

    def map_cos(self, expr):
        return expr.__class__(self.rec(expr.expr))