from vodes.symbolic.symbols import MachineError
from vodes.symbolic.interval import Interval
from pymbolic.mapper import RecursiveMapper

class ErrorMapper(RecursiveMapper):
    def __init__(self):
        if self is ErrorMapper:
            raise TypeError("ErrorAnalysis may not be directly instantiated")

    def map_call(self, expr):
        return self.rec(expr.function)(*[self.rec(par) for par in expr.parameters])
        
    # TODO : Constants can also introduce rounding error
    def map_constant(self, expr):
        return expr

class IntervalMapper(ErrorMapper):
    def map_variable(self, expr):
        return Interval(expr - expr * MachineError(), expr + expr * MachineError())

    def map_sum(self, expr):
        return sum(self.rec(child) for child in expr.children) * Interval(1 - MachineError(), 1 + MachineError())

    def map_product(self, expr):
        from pytools import product
        return product(self.rec(child) for child in expr.children) * Interval(1 - MachineError(), 1 + MachineError())

    def map_quotient(self, expr):
        return (self.rec(expr.numerator) / self.rec(expr.denominator)) * Interval(1 - MachineError(), 1 + MachineError())

    def map_power(self, expr):
        return (self.rec(expr.base) ** self.rec(expr.exponent)) * Interval(1 - MachineError(), 1 + MachineError())

    # TODO : Constants can also introduce rounding error
    def map_interval(self, expr):
        return expr


class AffineMapper(ErrorMapper):
    def map_variable(self, expr):
        pass

    def map_sum(self, expr):
        pass

    def map_product(self, expr):
        pass

    def map_quotient(self, expr):
        pass