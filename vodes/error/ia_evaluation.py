from pymbolic.mapper.evaluator import EvaluationMapper
from vodes.symbolic.symbols import Noise
from vodes.symbolic.interval import IntervalEvaluator
from vodes.symbolic.affine import Affine
from pymbolic.primitives import Expression, Sum, Product, Variable, is_constant
from pymbolic.mapper.stringifier import PREC_NONE, StringifyMapper

class MachineError(Variable):
    def __init__(self):
         super().__init__("e")
        
    def __lt__(self, other):
        if other >= 1:
            return True
        elif other <=0:
            return False
        else:
            return NotImplemented

class MachineEvaluator(IntervalEvaluator):
    def _eval_bound(self, expr):
        return MachineEvaluator(context=self.context, in_bounds=True)(expr)

    def _conversion(self, iv):
        # It is NOT possible to solve symbolic variables with interval functions -> convert to affine
        if is_constant(iv.low) and is_constant(iv.up):
            return iv

        #Affine arithmetic: concepts and applications
        return Affine(values=(
            (iv.up + iv.low)/2,
            (iv.up - iv.low)/2
        ), noises=(
             Noise(),
        ))

    def map_variable(self, expr):
        # we do not substitute machine error
        if not 'e' in self.context and isinstance(expr,MachineError):
            return expr
        else:
            return super().map_variable(expr)

    def map_interval(self, expr):
        res = super().map_interval(expr)
        return self._conversion(res)
        
    def _map_iv_product(self, expr):
        return self._conversion(super()._map_iv_product(expr))

    def _map_iv_quotient(self, expr):
        return self._conversion(super()._map_iv_quotient(expr))

    def _map_iv_sum(self, expr):
        return self._conversion(super()._map_iv_sum(expr))

class BoundEvaluator(EvaluationMapper):
    def map_variable(self, expr):
        # we do not substitute machine error
        if not 'e' in self.context and isinstance(expr,MachineError):
            return expr
        else:
            return super().map_variable(expr)
