from vodes.symbolic.expressions.primitives import Subtraction
from vodes.symbolic.expressions.trigonometric import sin, cos
from vodes.symbolic.expressions.nthroot import NthRoot
from vodes.symbolic.expressions.bounded import MachineError
from vodes.symbolic.expressions.interval import Interval

from pymbolic.mapper import RecursiveMapper
from pymbolic.primitives import Power, Variable

class ErrorMapper(RecursiveMapper):
    def __init__(self):
        if self is ErrorMapper:
            raise TypeError("ErrorAnalysis may not be directly instantiated")

    def map_call(self, expr):
        return self.rec(expr.function)(*[self.rec(par) for par in expr.parameters])
        
    # TODO : Constants can also introduce rounding error
    def map_constant(self, expr):
        return expr

    def map_rational(self, expr):
        return expr

class IntervalMapper(ErrorMapper):
    def map_variable(self, expr):
        return Interval(expr - expr * MachineError(), expr + expr * MachineError())

    def map_sum(self, expr):
        assert len(expr.children) == 2

        return sum(self.rec(child) for child in expr.children) * Interval(1 - MachineError(), 1 + MachineError())

    def map_sub(self, expr):
        assert len(expr.children) == 2

        return Subtraction(tuple([self.rec(child) for child in expr.children])) * Interval(1 - MachineError(), 1 + MachineError())

    def map_product(self, expr):
        from pytools import product

        assert len(expr.children) == 2

        return product(self.rec(child) for child in expr.children) * Interval(1 - MachineError(), 1 + MachineError())

    def map_quotient(self, expr):
        return (self.rec(expr.numerator) / self.rec(expr.denominator)) * Interval(1 - MachineError(), 1 + MachineError())

    def map_power(self, expr):
        return Power(self.rec(expr.base), self.rec(expr.exponent)) * Interval(1 - MachineError(), 1 + MachineError())

    # TODO : Constants can also introduce rounding error
    def map_interval(self, expr):
        return expr

    def map_nthroot(self, expr):
        return NthRoot(self.rec(expr.expr),expr.n) * Interval(1 - MachineError(), 1 + MachineError())

    def map_sin(self, expr):
        return sin(self.rec(expr.expr)) * Interval(1 - MachineError(), 1 + MachineError())

    def map_cos(self, expr):
        return cos(self.rec(expr.expr)) * Interval(1 - MachineError(), 1 + MachineError())

class AffineMapper(ErrorMapper):
    def __init__(self):
        super().__init__()

        self.__id = 0
        self.__noises = []
        self.__exact = {}

        self.__ids = {}

    @property
    def noises(self):
        """Get the noise variables that were created using this mapper"""
        return self.__noises

    @property
    def exact(self):
        """Get the substitutions of all created noise variables to zero"""
        return self.__exact

    @property
    def ids(self):
        """Get the mapping from ids to expressions"""
        return self.__ids
        
    def __error(self, expr):
        expr_id = str(expr)

        if not expr_id in self.__ids:
            # TODO : This would need to be tweaked in case of parallelization
            t = self.__id
            self.__id = self.__id + 1

            self.__ids[expr_id] = t

        eps = Variable(name=f'eps_{self.__ids[expr_id]}')

        if not eps.name in self.__exact:
            self.__noises.append(eps)
            self.__exact[eps.name] = 0

        return eps

    def map_variable(self, expr):
        return expr * (1 + self.__error(expr))

    def map_sum(self, expr):
        return sum(self.rec(child) for child in expr.children) * (1 + self.__error(expr))

    def map_sub(self, expr):
        assert len(expr.children) == 2

        return Subtraction(tuple([self.rec(child) for child in expr.children])) * (1 + self.__error(expr))

    def map_product(self, expr):
        from pytools import product
        return product(self.rec(child) for child in expr.children) * (1 + self.__error(expr))

    def map_quotient(self, expr):
        return (self.rec(expr.numerator) / self.rec(expr.denominator)) * (1 + self.__error(expr))

    def map_power(self, expr):
        return Power(self.rec(expr.base),self.rec(expr.exponent)) * (1 + self.__error(expr))

    # TODO : Constants can also introduce rounding error
    def map_interval(self, expr):
        return expr

    def map_nthroot(self, expr):
        return NthRoot(self.rec(expr.expr),expr.n) * (1 + self.__error(expr))

    def map_sin(self, expr):
        return sin(self.rec(expr.expr)) * (1 + self.__error(expr))

    def map_cos(self, expr):
        return cos(self.rec(expr.expr)) * (1 + self.__error(expr))