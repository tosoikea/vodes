
from pymbolic.mapper import RecursiveMapper
from pymbolic.primitives import Variable

import sympy

class ExactPymbolicToSympyMapper(RecursiveMapper):
    """"""
    
    @property
    def sym(self):
        return sympy
    
    def map_constant(self, expr):
        return self.sym.sympify(expr)

    def map_rational(self,expr):
        return self.sym.Rational(expr.numerator, expr.denominator)

    def map_variable(self, expr):
        return self.sym.Symbol(expr.name)
        
    def map_subscript(self, expr):
        return self.sym.Indexed(
            self.sym.IndexedBase(self.rec(expr.aggregate)),
            *tuple(self.rec(i) for i in expr.index_tuple)
            )

    def map_call(self, expr):
        if isinstance(expr.function, Variable):
            func_name = expr.function.name
            try:
                func = getattr(self.sym.functions, func_name)
            except AttributeError:
                func = self.sym.Function(func_name)
            return func(*[self.rec(par) for par in expr.parameters])
        else:
            self.raise_conversion_error(expr)

    def map_sum(self, expr):
        return self.sym.Add(
            *[self.rec(child) for child in expr.children]
        ) 

    def map_product(self, expr):
        return self.sym.Mul(
            *[self.rec(child) for child in expr.children]
        ) 
        
    def map_quotient(self, expr):
        return self.sym.Mul(
            self.rec(expr.numerator),
            self.sym.Pow(
                self.rec(expr.denominator),
                -1
            )
        )

    def map_power(self, expr):
        return self.sym.Pow(
            self.rec(expr.base),
            self.rec(expr.exponent)
        )  


    def map_derivative(self, expr):
        return self.sym.Derivative(self.rec(expr.child),
                *[self.sym.Symbol(v) for v in expr.variables])