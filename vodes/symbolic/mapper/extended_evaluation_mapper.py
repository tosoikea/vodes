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
    def map_sub(self,expr):
        from functools import reduce

        return reduce(
            lambda r,x: r - x,
            [ self.rec(child) for child in expr.children ] 
        )

    def map_interval(self,expr):
        from vodes.symbolic.expressions.interval import Interval
        
        return Interval(
            self.rec(expr.low),
            self.rec(expr.up)
        )
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
    