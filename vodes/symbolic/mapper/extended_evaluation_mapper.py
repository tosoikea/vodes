from vodes.symbolic.expressions.nthroot import NthRoot
from pymbolic.mapper.evaluator import EvaluationMapper
from math import pi, sin, cos

def evaluate(expression, context=None, float:bool=False):
    if context is None:
        context = {}
    res = ExtendedEvaluationMapper(context)(expression)

    # Two iterations of solver, if symbolic values are used for evaluation.
    # This allows to push the floating calculations further up.
    if float:
        return evaluate(res, context=context)
    else:
        return res

class ExtendedEvaluationMapper(EvaluationMapper):
    ###
    # FUNCTIONS
    ###
    def map_nthroot(self,expr):
        return self.rec(self.rec(expr.expr)**(1/expr.n))

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
    