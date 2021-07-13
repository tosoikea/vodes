from sys import intern
from vodes.symbolic.symbols import Boundary, BoundedVariable
from pymbolic.mapper import RecursiveMapper
from pymbolic.mapper.stringifier import StringifyMapper
from pymbolic.primitives import Expression, Product, Sum, Variable


class Abs(Expression):
    """Class to represent the absolute function as an expression.

    Args:
        expression: The expression to be wrapped within the absolute function.

    Attributes:
        expression: The wrapped expression.
    """
    init_arg_names = ("expression",)

    def __init__(self, expression):
        assert(not expression is None)

        self.expression = expression

    def __getinitargs__(self):
        return (self.expr,)

    @property
    def expr(self):
        return self.expression

    def make_stringifier(self, originating_stringifier=None):
        return AbsoluteStringifyMapper()

    mapper_method = intern("map_absolute")

class AbsoluteStringifyMapper(StringifyMapper):
    def map_absolute(self, expr, enclosing_prec, *args, **kwargs):
        return f'|{expr.expr}|'

class AbsoluteEvaluator(RecursiveMapper):
    """Class to evaluate expressions on the basis of the absolute function. It supports free variables on the basis of provided bounds.

    Args:
        absolute (bool): Is the evaluator currently within an abs(...) call? If yes, all further subexpressions are evaluated as absolutes.
        b (Boundary): Specify a boundary for the expression.
    """
    def __init__(self, absolute=False, b:Boundary=None):
        self.absolute = absolute
        self.boundary = b

    def map_absolute(self, expr:Abs) -> Expression:
        return AbsoluteEvaluator(absolute=True)(expr.expr)

    def map_constant(self, expr) -> Expression:
        return abs(expr)

    def map_variable(self, expr:Variable) -> Expression:
        boundaries = []

        # TODO : Assert correctness
        if isinstance(expr, BoundedVariable):
            boundaries.append(expr.bound)

        if not (self.boundary is None):
            boundaries.append(self.boundary)

        for b in boundaries:
            # [0,x] => positive
            if b.lower.value >= 0:
                return expr
            # [y,0] => negative
            elif b.upper.value <= 0:
                return (-1) * expr
            
        raise ValueError(f"Cannot evaluate symbol {expr}")