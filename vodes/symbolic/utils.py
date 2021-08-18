from typing import List, Tuple
from vodes.symbolic.expressions.infinity import Infinity, NegativeInfinity
from pymbolic.primitives import Expression

_prio_types = (Infinity, NegativeInfinity)

def gt(o,v) -> bool:
    if isinstance(o, _prio_types):
        return o > v
    
    if isinstance(v, _prio_types):
        return v < o

    if isinstance(o,Expression):
        from pymbolic import evaluate
        return evaluate(o.gt(v))

    if isinstance(v,Expression):
        from pymbolic import evaluate
        return evaluate(v.lt(o))
    
    return o > v

def ge(o,v) -> bool:
    if isinstance(o, _prio_types):
        return o >= v
    
    if isinstance(v, _prio_types):
        return v <= o

    if isinstance(o,Expression):
        from pymbolic import evaluate
        return evaluate(o.ge(v))

    if isinstance(v,Expression):
        from pymbolic import evaluate
        return evaluate(v.le(o))
    
    return o >= v

def lt(o,v) -> bool:
    return gt(v,o)

def le(o,v) -> bool:
    return ge(v,o)

def eq(o,v) -> bool:
    if isinstance(o, _prio_types) or isinstance(v, _prio_types):
        return o == v

    if isinstance(o,Expression):
        from pymbolic import evaluate
        return evaluate(o.eq(v))

    if isinstance(v,Expression):
        from pymbolic import evaluate
        return evaluate(v.eq(o))
    
    return o == v

def minimum(os:list):
    if len(os) == 0:
        raise ValueError("Cannot determine minimum for empty list.")

    res = os[0]
    for i in range(len(os) - 1):
        if lt(os[i+1],res):
            res = os[i+1]

    return res
    
def maximum(os:list):
    if len(os) == 0:
        raise ValueError("Cannot determine maximum for empty list.")

    res = os[0]
    for i in range(len(os) - 1):
        if gt(os[i+1],res):
            res = os[i+1]

    return res


def merge(ls:list,rs:list) -> List[Tuple[Expression,Expression,object]]:
    """Merges two list of bounded expressions to tuples made up of a pair of items from the lists and the shared domain."""
    res = []

    for l in ls:
        for r in rs:
            combined = l.bound.intersect(r.bound)

            if combined is None:
                continue

            res.append(
                (
                    l.expr,
                    r.expr,
                    combined
                )
            )

    return res


