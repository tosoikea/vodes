from pymbolic.mapper import RecursiveMapper
from pymbolic.primitives import Sum, Product


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
    def __split(self,f, expr):   
        if len(expr.children) <= 2:
            return expr
            
        children = [
            expr.children[0],
            expr.children[1]
        ]

        for i in range (2, len(expr.children)):
            children[0] = f(children[0],children[1])
            children[1] = expr.children[i]

        return f(*children)

    def map_sum(self, expr): 
        return self.__split(lambda a,b : Sum((a,b)), expr)

    def map_product(self, expr):
        return self.__split(lambda a,b : Product((a,b)), expr)

    # always binary (numerator, denominator)
    def map_quotient(self, expr):
        return expr

    def map_constant(self, expr):
        return expr
        
    def map_variable(self, expr):
        return expr