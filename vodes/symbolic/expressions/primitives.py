from sys import intern

# Expression Library
from pymbolic.primitives import Expression, Variable, is_constant, is_nonzero, is_valid_operand

# Expression Mapper
from pymbolic.mapper import RecursiveMapper
from pymbolic.mapper.stringifier import StringifyMapper, PREC_NONE

class ExtendedExpression(Expression):
    def __sub__(self, other):
        if not is_valid_operand(other):
            return NotImplemented

        if is_nonzero(other):
            if self:
                if isinstance(other, Subtraction):
                    return Subtraction((self,) + other.children)
                else:
                    return Subtraction((self, other))
            else:
                return other
        else:
            return self

    def __rsub__(self, other):
        if not is_constant(other):
            return NotImplemented

        if is_nonzero(other):
            return Subtraction((other,self))
        else:
            return -self

class Subtraction(ExtendedExpression):
    """Class to represent subtraction. In Pymbolic this is handled using an addition and multiplication with -1. However, this would introduce further errors into the error analysis.

    Args:
        children: The children to subtract from each other.
    """
    init_arg_names = ("children",)

    def __init__(self, children):
        assert isinstance(children, tuple)

        self.children = children

    def __getinitargs__(self):
        return (self.children,)

    @property
    def args(self):
        return self.children

    def __sub__(self, other):
        if not is_valid_operand(other):
            return NotImplemented

        if isinstance(other, Subtraction):
            return Subtraction(self.children + other.children)
        if not other:
            return self
        return Subtraction(self.children + (other,))

    def __rsub__(self, other):
        if not is_constant(other):
            return NotImplemented

        if isinstance(other, Subtraction):
            return Subtraction(other.children + self.children)
        if not other:
            return self
        return Subtraction((other,) + self.children)

    def make_stringifier(self, originating_stringifier=None):
        return PrimitiveStringifyMapper()

    mapper_method = intern("map_sub")

class PrimitiveStringifyMapper(StringifyMapper):
    def map_sub(self, expr, enclosing_prec, *args, **kwargs):
        return self.parenthesize_if_needed(
                self.join_rec(" - ", expr.children, PREC_NONE, *args, **kwargs),
                enclosing_prec, PREC_NONE)