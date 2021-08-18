from sys import intern
from typing import List

# Custom Expression library
from vodes.symbolic.utils import lt,gt,eq,le,ge

# Expression Library
from pymbolic.primitives import Expression, Variable, Power

# Expression Mapper
from pymbolic.mapper.stringifier import StringifyMapper


class Domain:
    """Class to represent the domain of an expression.
    
    Args:
        start: Beginning of the domain.
        end: End of the domain.
        left_open: If this value is true, the start is not part of the domain.
        right_open: If this value is true, the end is not part of the domain."""

    def __init__(self, start, end=None, left_open:bool=False, right_open:bool=False):
        assert not (start is None)

        if end is None:
            end = start

        assert le(start,end)
        
        if eq(start,end) and right_open:
            raise ValueError(f"{start} equals end and right side is open. This is invalid.")

        self.start = start
        self.end = end
        self.left_open=left_open
        self.right_open=right_open

    ###
    # Domain functions
    ###
    def contains(self, item):
        """Determines if the supplied item is part of the domain."""
        res = False

        if self.left_open:
            res = lt(self.start,item)
        else:
            res = le(self.start,item)

        if self.right_open:
            res = gt(res and self.end,item)
        else:
            res = ge(res and self.end,item)

        return res

    def intersect(self, other):
        """Constructs the intersection of this instance with another domain."""

        if not isinstance(other,Domain):
            raise TypeError(f"Invalid type {type(other)} for intersection of domains.")

        lv = self.start
        lo = self.left_open

        if gt(other.start, lv):
            lv = other.start
            lo = other.left_open
        elif eq(other.start, lv):
            lo = other.left_open or self.left_open

        rv = self.end
        ro = self.right_open

        if lt(other.end, rv):
            rv = other.end
            ro = other.right_open
        elif eq(other.end, rv):
            ro = other.right_open and self.right_open

        if lt(rv, lv):
            return None
        elif eq(rv, lv) and (ro or lo):
            return None
        else:
            return Domain(lv, rv, lo, ro)

    def union(self, other):
        """Constructs the union of this instance with another domain."""

        if not isinstance(other,Domain):
            raise TypeError(f"Invalid type {type(other)} for union of domains.")

        lv = self.start
        lo = self.left_open
        
        if lt(other.start, lv):
            lv = other.start
            lo = other.left_open
        elif eq(other.start, lv):
            lo = other.left_open and self.left_open

        rv = self.end
        ro = self.right_open

        if gt(other.end, rv):
            rv = other.end
            ro = other.right_open
        elif eq(other.end, rv):
            ro = other.right_open and self.right_open
        
        return Domain(lv, rv, lo, ro)

    def difference(self, other):
        """Constructs the difference (set minus) of this instance with a list of other domains."""
        pass

    ####
    # Built-In Functions
    ####   
    def __eq__(self, o: object) -> bool:
        if isinstance(o, Domain):
            return (
                o.start == self.start and o.left_open == self.left_open and
                o.end == self.end and o.right_open == self.right_open
            )

        return False

    def __hash__(self):
        return hash((self.start, self.end, self.left_open, self.right_open))

    def __str__(self) -> str:
        left = f'({self.start}' if self.left_open else f'[{self.start}'
        right = f'{self.end})' if self.right_open else f'{self.end}]'

        return f'{left},{right}'

    def __repr__(self) -> str:
        return str(self)

class BoundedExpression(Expression):
    init_arg_names = ("boundary","expression",)

    def __init__(self,expression,boundary:Domain):
        self.boundary = boundary
        self.expression = expression

    def __getinitargs__(self):
        return self.expression, self.boundary

    @property
    def expr(self):
        return self.expression

    @property
    def bound(self):
        return self.boundary

    def make_stringifier(self, originating_stringifier=None):
        return BoundedExpressionStringifyMapper()

    mapper_method = intern("map_bounded_expression")

class BoundedExpressionStringifyMapper(StringifyMapper):
    def map_bounded_expression(self, expr, enclosing_prec, *args, **kwargs):
        return f'B({expr.bound},{expr.expr})'

class BoundedVariable(Variable):
    init_arg_names = ("boundary","name",)

    def __init__(self,name:str,boundary:Domain):
        super().__init__(name=name)
        self.__boundary = boundary

    def __getinitargs__(self):
        return (self.bound, self.name, )

    @property
    def bound(self):
        return self.__boundary

class MachineError(BoundedVariable):
    def __init__(self,min_precision:int=1,max_precision:int=None):
        assert (min_precision > 0)
        max_eps = Power(base=2, exponent=-min_precision)
        
        min_eps = 0
        if not max_precision is None:
            assert (max_precision >= min_precision)
            min_eps = Power(base=2, exponent=-max_precision)

        super().__init__(
            name="eps",
            boundary=Domain(
                # open if value unreachable
                start=min_eps,
                left_open=min_eps == 0,
                # open if value unreachable
                end=max_eps,
                right_open=max_eps == 1
            )
        )

class Noise(BoundedVariable):
    def __init__(self,index:str):
        super().__init__(
            name=f'eps_{index}',
            boundary=Domain(
                # open if value unreachable
                start=-MachineError(),
                # open if value unreachable
                end=MachineError()
            )
        )

class DummyVariable(BoundedVariable):
    def __init__(self):
        super().__init__(
            name="dummy",
            boundary=Domain(
                start=-1,
                left_open=True,
                end=1,
                right_open=True
            )
        )