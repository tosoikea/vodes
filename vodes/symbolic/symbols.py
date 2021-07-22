from operator import attrgetter
from typing import List
from pymbolic.primitives import Expression, Variable
from sys import intern
from pymbolic.mapper.stringifier import PREC_NONE, StringifyMapper
from vodes.symbolic.power import Power

class BoundedValue:
    def __init__(self, value, open:bool):
        assert not (value is None)
        
        self.value = value
        self.open = open

    def __eq__(self, o: object) -> bool:
        if isinstance(o, BoundedValue):
            return o.value == self.value and o.open == self.open

        return False

    def __hash__(self):
        return hash((self.value, self.open))

    def __str__(self) -> str:
        if self.open:
            return f'B({self.value})'
        else:
            return f'B[{self.value}]'

    def __repr__(self) -> str:
        return str(self)

class Boundary:
    def __init__(self, lower:BoundedValue, upper:BoundedValue):
        assert lower
        assert upper

        self.lower = lower
        self.upper = upper

    def __eq__(self, o: object) -> bool:
        if isinstance(o, Boundary):
            return self.lower == o.lower and self.upper == o.upper

        return False

    def __hash__(self):
        return hash((self.lower, self.upper))

    def contains(self, item):
        res = False

        if self.lower.open:
            res = self.lower.value < item
        else:
            res = self.lower.value <= item

        if self.upper.open:
            res = res and self.upper.value > item
        else:
            res = res and self.upper.value >= item

        return res

    def __str__(self) -> str:
        return f'{"(" if self.lower.open else "["}{self.lower.value},{self.upper.value}{")" if self.upper.open else "]"}'

    def __repr__(self) -> str:
        return str(self)

    def intersect(self, other):
        if not isinstance(other,Boundary):
            raise TypeError(f"Invalid type {type(other)} for intersection of boundary intervals.")

        l = self.lower
        if other.lower.value > l.value:
            l = other.lower
        elif other.lower.value == l.value:
            l = BoundedValue(value = l.value, open=l.open or other.lower.open)

        r = self.upper
        if other.upper.value < r.value:
            r = other.upper
        elif other.upper.value == r.value:
            r = BoundedValue(value = r.value, open=r.open or other.upper.open)

        if r.value < l.value:
            return None
        elif r.value == l.value and (r.open or l.open):
            return None
        else:
            return Boundary(lower=l,upper=r)

    def union(self, other):
        if not isinstance(other,Boundary):
            raise TypeError(f"Invalid type {type(other)} for union of boundary intervals.")

        l = self.lower
        if other.lower.value < l.value:
            l = other.lower
        elif other.lower.value == l.value:
            l = BoundedValue(value = l.value, open=l.open and other.lower.open)

        r = self.upper
        if other.upper.value > r.value:
            r = other.upper
        elif other.upper.value == r.value:
            r = BoundedValue(value = r.value, open=r.open and other.upper.open)
        
        return Boundary(lower=l,upper=r)


    # TODO : Expects other sorted by lower boundary
    def difference(self, other):
        if not isinstance(other,list):
            raise TypeError(f"Invalid type {type(other)} for difference of boundary intervals.")

        other = sorted(other, key=attrgetter("lower"))

        res = []

        l = self.lower

        for b in other:
            if l < b.lower:
                res.append(
                    Boundary(
                        lower=l,
                        upper=BoundedValue(
                            value=b.lower.value,
                            open=not b.lower.open
                        )
                    )
                )
            l = b.upper

        return res


class BoundedExpression(Expression):
    init_arg_names = ("boundary","expression",)

    def __init__(self,expression,boundary:Boundary):
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

    def __init__(self,name:str,boundary:Boundary):
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
            boundary=Boundary(
                # open if value unreachable
                lower=BoundedValue(value=min_eps,open=min_eps == 0),
                # open if value unreachable
                upper=BoundedValue(value=max_eps,open=max_eps == 1)
            )
        )

class DummyVariable(BoundedVariable):
    def __init__(self):
        super().__init__(
            name="dummy",
            boundary=Boundary(
                lower=BoundedValue(value=-1,open=True),
                upper=BoundedValue(value=1,open=True)
            )
        )

# TODO : Move to BoundedVariable
class Noise(Variable):
    glb_index = 0

    def __init__(self,index=None):
        # TODO : Find better way for unique index -> overflow, race condition etc.
        if index is None:
            index = Noise.glb_index
            Noise.glb_index += 1

        self.index = index
        super().__init__(f"e_{index}")
        
    def __lt__(self, other):
        if other > 1:
            return True
        elif other <-1:
            return False
        else:
            return NotImplemented