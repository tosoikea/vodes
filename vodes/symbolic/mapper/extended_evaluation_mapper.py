from vodes.symbolic.expressions.nthroot import NthRoot
from pymbolic.mapper.evaluator import EvaluationMapper
from math import pi, sin, cos

def evaluate(expression, context=None):
    if context is None:
        context = {}
    return ExtendedEvaluationMapper(context)(expression)

class ExtendedEvaluationMapper(EvaluationMapper):
    ###
    # FUNCTIONS
    ###
    def map_nthroot(self,expr):
        return evaluate(self.rec(expr.expr)**(1/expr.n))

    def map_sin(self, expr):
        return sin(
            self.rec(expr.expr)
        )

    def map_cos(self, expr):
        return cos(
            self.rec(expr.expr)
        )
    ###
    # CONSTANTS
    ###
    def map_pi(self,expr):
        return pi
    