
from typing import final
from vodes.symbolic.expressions.interval import Interval
from pymbolic.mapper import RecursiveMapper
from pymbolic.primitives import Variable, Quotient

import sympy

class ExactSympyToPymbolicMapper:
    def __init__(self, context:dict=None) -> None:
        super().__init__()
        self._context = {} if context is None else context

    def __call__(self, expr):
        return self.__after_walk(
            expr=expr,
            f=self.__to_pymbolic
        )

    def __after_walk(self,expr,f):
        from sympy.core.basic import Basic
        from sympy.sets.setexpr import SetExpr
        from sympy.sets.sets import FiniteSet, Interval as SymInterval

        if not isinstance(expr,Basic):
            return expr
            
        if isinstance(expr, SetExpr):
            if isinstance(expr.set,FiniteSet):
                return Interval(
                    self.__after_walk(expr.set.args[0], f)
                )
            elif isinstance(expr.set,SymInterval):
                return Interval(
                    self.__after_walk(expr.set.start,f),
                    self.__after_walk(expr.set.end,f)
                )
            else:
                raise ValueError("{0}.{1}".format(expr.set.__class__.__module__,expr.set.__class__.__name__))

        children = []

        for arg in expr.args:
            children.append(self.__after_walk(arg,f))
        
        return f(expr,children)

    def __to_pymbolic(self,expr,children):
        from pymbolic.primitives import Product, Sum, Variable, Power
        from vodes.symbolic.expressions.constants import Pi
        from vodes.symbolic.expressions.infinity import Infinity, NegativeInfinity
        from vodes.symbolic.expressions.nthroot import NthRoot
        from vodes.symbolic.expressions.absolute import Abs
        from sympy.core.numbers import Pi as SymPi
        from sympy.functions.elementary.complexes import Abs as SymAbs


        ## (0) Variables
        if expr.is_symbol:
            name = str(expr)

            if name in self._context:
                return self._context[name]
            else:
                return Variable(str(expr))

        ## (1) Constants

        # Integer, Rational, Zero, One, NegativeOne, Half
        elif expr.is_Rational:
            if expr.p == 0:
                return 0
            elif expr.q == 1:
                return expr.p
            else:
                return Quotient(expr.p, expr.q)
        elif expr.is_infinite:
            if expr.is_extended_positive:
                return Infinity()
            else:
                return NegativeInfinity()
        elif isinstance(expr,SymPi):
            return Pi()

        ## (2) Unary 

        ## (3) Arithmetic
        elif expr.is_Mul:
            return Product(tuple(children))
        elif expr.is_Add:
            return Sum(tuple(children))
        elif expr.is_Pow:
            # x^(0.5) => sqrt(x), x^(1.5) => sqrt(x)^3
            base = children[0]
            exponent = children[1]

            if isinstance(exponent,Quotient) and exponent.den > 1:
                return NthRoot(expression=base, n=exponent.den) ** exponent.num
            
            return Power(base,exponent)
        elif isinstance(expr, SymAbs):
            return Abs(children[0])
        else:
            raise ValueError("{0}.{1}".format(expr.__class__.__module__,expr.__class__.__name__))


class ExactPymbolicToSympyMapper(RecursiveMapper):
    """Mapper for converting pymbolic expression to exact representations within sympy (without evaluation)"""
    def __init__(self, context:dict=None) -> None:
        super().__init__()
        self._context = {} if context is None else context

    def map_foreign(self, expr, *args, **kwargs):
        """Mapper method dispatch for non-:mod:`pymbolic` objects."""
        try:
            return super().map_foreign(expr,*args,**kwargs)
        except ValueError:
            raise ValueError(f"{self.__class__} encountered invalid foreign object: {repr(expr)}({expr.__class__})")

    @property
    def sym(self):
        return sympy
    
    def map_constant(self, expr):
        return self.sym.sympify(expr)

    def map_rational(self,expr):
        return self.map_quotient(expr)

    def map_interval(self,expr):
        from sympy.sets.setexpr import SetExpr

        return SetExpr(
            self.sym.sets.Interval(
                self.rec(expr.low),
                self.rec(expr.up)
            )
        )

    def map_variable(self, expr):
        if expr.name in self._context:
            return self._context[expr.name]
        else:
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

    def map_sub(self, expr):
        children = []
        for i in range(len(expr.children)):
            child = expr.children[i]
            children.append(
                self.rec(child) if i == 0 else self.sym.Mul(-1,self.rec(child))
            )

        return self.sym.Add(*children)
        
    def map_product(self, expr):
        return self.sym.Mul(
            *[self.rec(child) for child in expr.children]
        ) 
        
    def map_quotient(self, expr):
        return self.sym.Mul(
                self.rec(expr.num),
                self.sym.Pow(
                    self.rec(expr.den),
                    -1
                )
            )

    def map_power(self, expr):
        return self.sym.Pow(
            self.rec(expr.base),
            self.rec(expr.exponent),
            evaluate=False
        )  


    def map_derivative(self, expr):
        return self.sym.Derivative(self.rec(expr.child),
                *[self.sym.Symbol(v) for v in expr.variables])

    ###
    # FUNCTIONS
    ###
    def map_nthroot(self,expr):
        return self.sym.Pow(
            self.rec(expr.expr),
            self.map_quotient(Quotient(1,expr.n))
        )

    def map_sin(self, expr):
        return self.sym.sin(
            self.rec(expr.expr)
        )

    def map_cos(self, expr):
        return self.sym.cos(
            self.rec(expr.expr)
        )

    def map_absolute(self, expr):
        return self.sym.functions.elementary.complexes.Abs(
            self.rec(expr.expr)
        )

    ###
    # CONSTANTS
    ###
    def map_pi(self,expr):
        return self.sym.pi

class ExactPymbolicToMathMapper(RecursiveMapper):
    """Mapper for converting pymbolic expression to representations within mpmath (with evaluation)"""
    def __init__(self, precision:int, context:dict=None) -> None:
        super().__init__()
        self._prec = precision
        self._context = {} if context is None else context

    def __call__(self, expr, *args, **kwargs):
        from pymbolic.primitives import Expression
        from mpmath import mp

        try:
            method = getattr(self, expr.mapper_method)
        except AttributeError:
            if isinstance(expr, Expression):
                return self.handle_unsupported_expression(
                        expr, *args, **kwargs)
            else:
                return self.map_foreign(expr, *args, **kwargs)

        t = mp.prec
        try:
            mp.prec = self._prec
            return method(expr, *args, **kwargs)
        finally:
            mp.prec = t

    def map_foreign(self, expr, *args, **kwargs):
        """Mapper method dispatch for non-:mod:`pymbolic` objects."""
        try:
            return super().map_foreign(expr,*args,**kwargs)
        except ValueError:
            raise ValueError(f"{self.__class__} encountered invalid foreign object: {repr(expr)}({expr.__class__})")
    
    def map_constant(self, expr):
        from mpmath import mpf
        return mpf(expr)

    def map_rational(self,expr):
        return self.map_quotient(expr)

    def map_variable(self, expr):
        from pymbolic.mapper.evaluator import UnknownVariableError
        from mpmath import mpf

        try:
            subs = self._context[expr.name]
        except KeyError:
            raise UnknownVariableError(expr.name)
        
        return mpf(subs)

    def map_sum(self, expr):
        return sum([self.rec(child) for child in expr.children])

    def map_sub(self, expr):
        children = []
        for i in range(len(expr.children)):
            child = expr.children[i]
            children.append(
                self.rec(child) if i == 0 else (-1) * self.rec(child)
            )

        return sum(children)
        
    def map_product(self, expr):
        from pytools import product
        return product([self.rec(child) for child in expr.children])
        
    def map_quotient(self, expr):
        return self.rec(expr.num) / self.rec(expr.den)

    def map_power(self, expr):
        return self.rec(expr.base) ** self.rec(expr.exponent)

    ###
    # FUNCTIONS
    ###
    def map_nthroot(self,expr):
        from mpmath import nthroot
        return nthroot(self.rec(expr),expr.n)

    def map_sin(self, expr):
        from mpmath import sin
        return sin(expr.expr)

    def map_cos(self, expr):
        from mpmath import cos
        return cos(expr.expr)

    ###
    # CONSTANTS
    ###
    def map_pi(self,expr):
        from mpmath import pi
        return pi
