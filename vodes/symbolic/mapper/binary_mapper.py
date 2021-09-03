import logging
from vodes.symbolic.expressions.nthroot import NthRoot
from vodes.symbolic.expressions.primitives import Subtraction

from pymbolic.mapper import RecursiveMapper
from pymbolic.primitives import Sum, Product, Quotient, Power


class BinaryMapper(RecursiveMapper):
    """Traverses an expression to convert it into a binary tree. 
    This simplifies e.g. later analysis, as every addition adds an error.
    
    Example usage:
    .. doctest::
        >>> from pymbolic import var
        >>> x = var("x")
        >>> y = var("y")
        >>> u = 5 + x + y
        >>> repr(u)
        "Sum((5, Variable('x'), Variable('y')))"
        >>> from vodes.symbolic.binary_mapper import BinaryMapper as BM
        >>> repr(BM()(u))
        "Sum((Sum((5, Variable('x'))), Variable('y')))"
    """  
    def __init__(self):
        self._logger = logging.getLogger(__name__)
    
    def __split(self,f, expr):   
        binary_children = [self.rec(child) for child in expr.children]

        if len(binary_children) <= 2:
            return f(*binary_children)

        children = [
            binary_children[0],
            binary_children[1]
        ]

        for i in range (2, len(binary_children)):
            children[0] = f(children[0],children[1])
            children[1] = binary_children[i]

        return f(*children)

    def map_sum(self, expr): 
        res =  self.__split(lambda a,b : Sum((a,b)), expr)
        self._logger.debug(f'{expr} -> {repr(res)}')
        return res

    def map_sub(self, expr): 
        res =  self.__split(lambda a,b : Subtraction((a,b)), expr)
        self._logger.debug(f'{expr} -> {repr(res)}')
        return res

    def map_product(self, expr):
        res = self.__split(lambda a,b : Product((a,b)), expr)
        self._logger.debug(f'{expr} -> {repr(res)}')
        return res

    # always binary (numerator, denominator)
    def map_rational(self, expr):
        return self.map_quotient(expr)

    # always binary (numerator, denominator)
    def map_quotient(self, expr):
        res = Quotient(
            self.rec(expr.num),
            self.rec(expr.den)
        )
        self._logger.debug(f'{expr} -> {repr(res)}')
        return res

    # always binary (base, exponent)
    def map_power(self, expr):
        return Power(
            self.rec(expr.base),
            self.rec(expr.exponent)
        )

    ## always binary (expr,n)
    def map_nthroot(self, expr):
        return NthRoot(
            self.rec(expr.expr),
            expr.n
        )
        
    # always unary
    def map_constant(self, expr):
        return expr
        
    # always unary
    def map_variable(self, expr):
        return expr
        
    # always unary
    def map_sin(self, expr):
        raise NotImplementedError()

    # always unary
    def map_cos(self, expr):
        raise NotImplementedError()
