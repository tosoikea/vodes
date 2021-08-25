from sys import intern
from typing import List

# Custom Expression library
from vodes.symbolic.utils import lt,gt,eq,le,ge
from vodes.symbolic.expressions.primitives import ExtendedExpression

# Expression Library
from pymbolic.primitives import Variable, Power

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

    def intersect(self, other) -> object:
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

    def union(self, other) -> list:
        """Constructs the union of this instance with another domain."""

        if not isinstance(other,Domain):
            raise TypeError(f"Invalid type {type(other)} for intersection of domains.")

        # (1) Is outside
        if self.start > other.end or self.end < other.start:
            return [self,other]
        # (2) Expand
        else: 
            (lv,lo) = (self.start,self.left_open)

            if lt(other.start,self.start):
                (lv,lo) = (other.start,other.left_open)
            elif eq(other.start,self.start):
                lo = other.left_open and self.left_open

            (rv,ro) = (self.end,self.right_open)

            if gt(other.end,self.end):
                (rv,ro) = (other.end,other.right_open)
            elif eq(other.end,self.end):
                lo = other.right_open and self.right_open

            return [Domain(lv,rv,lo,ro)]

    def difference(self, others:list) -> list:
        """Constructs the difference (set minus) of this instance with a list of other domains."""
        res = [self]

        for other in others:
            assert isinstance(other,Domain)
            
            tres = []
            for r in res:
                # (1) Is outside, no effect
                # a) ..] ... [..
                # b) )[.., ](.., )(..
                # c) ..](, ..)[, ..)(
                if is_outside(r,other):
                    tres.append(r)
                # (2) Splitting into two
                # This is given, when r is fully including other and contains more items
                elif (lt(r.start,other.start) or (eq(r.start,other.start) and other.left_open and not r.left_open)) and (gt(r.end,other.end) or (eq(r.end,other.end) and other.right_open and not r.right_open)):
                    # 1.
                    (rv,ro) = (other.start, not other.left_open)
                    tres.append(
                        Domain(
                            r.start,rv,r.left_open,ro
                        )
                    )

                    # 2.
                    (lv,lo) = (other.end, not other.right_open)
                    tres.append(
                        Domain(
                            lv,r.end,lo,r.right_open
                        )
                    )
                # (3) Splitting into single
                else:
                    lv,lo,rv,ro = None,None,None,None

                    if gt(other.start,r.start) or (eq(other.start,r.start) and other.left_open and not r.left_open):
                        (lv,lo) = (r.start,r.left_open)
                        rv = other.start
                        ro = not other.left_open
                    elif lt(other.end,r.end) or (eq(other.end,r.end) and other.right_open and not r.right_open):
                        (rv,ro) = (r.end,r.right_open)
                        lv = other.end
                        lo = not other.right_open
                    else:
                        continue
                        
                    tres.append(
                        Domain(
                            lv,rv,lo,ro
                        )
                    )

            res = tres

        return res

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

class BoundedExpression(ExtendedExpression):
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

####
# Utility Functions
####   
def is_outside(r:Domain,other:Domain) -> bool:
    return (gt(r.start,other.end) or lt(r.end,other.start)) or (eq(r.start,other.end) and (r.left_open or other.right_open)) or (eq(r.end,other.start) and (r.right_open or r.left_open))
    
def intersect(ls:list,rs:list):
    res = []

    if ls is None:
        ls = []
    
    if rs is None:
        rs = []

    for l in ls:
        for r in rs:
            combined = l.intersect(r)

            if combined is None:
                continue
        
            res.append(combined)
            
    return res 

