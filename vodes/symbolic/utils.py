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
        from vodes.symbolic.mapper.extended_evaluation_mapper import evaluate
        return evaluate(o.gt(v))

    if isinstance(v,Expression):
        from vodes.symbolic.mapper.extended_evaluation_mapper import evaluate
        return evaluate(v.lt(o))
    
    return o > v

def ge(o,v) -> bool:
    if isinstance(o, _prio_types):
        return o >= v
    
    if isinstance(v, _prio_types):
        return v <= o

    if isinstance(o,Expression):
        from vodes.symbolic.mapper.extended_evaluation_mapper import evaluate
        return evaluate(o.ge(v))

    if isinstance(v,Expression):
        from vodes.symbolic.mapper.extended_evaluation_mapper import evaluate
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
        from vodes.symbolic.mapper.extended_evaluation_mapper import evaluate
        return evaluate(o.eq(v))

    if isinstance(v,Expression):
        from vodes.symbolic.mapper.extended_evaluation_mapper import evaluate
        return evaluate(v.eq(o))
    
    return o == v

def compare(item1,item2) -> int:
    """Compare function for sorting (constant) boundaries"""
    if lt(item1.start,item2.start):
        return -1
    elif gt(item1.start,item2.start):
        return 1
    else:
        if lt(item1.end,item2.end):
            return -1
        elif gt(item1.end,item2.end):
            return 1
        else:
            return 0

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


def sort_bounded(ls:list):
    """Sort a list of bounded expressions"""
    from functools import cmp_to_key
    list.sort(ls,key=cmp_to_key(lambda item1,item2: compare(item1.bound,item2.bound)))
    

def merge(ls:list,rs:list) -> List[Tuple[Expression,Expression,object]]:
    """Merges two list of bounded expressions to tuples made up of a pair of items from the lists and the shared domain."""

    res = []

    sort_bounded(ls)
    sort_bounded(rs)
  
    for l in ls:
        for r in rs:
            combined = l.bound.intersect(r.bound)

            if combined is None:
                continue

            # Multiple domains contain the same item => Ambiguity
            if len(res) > 0 and res[-1][2].right_open == combined.left_open and eq(combined.start,res[-1][2].end):
                    # a) Single item domain
                    # Current one? => Ignore, information is already present within given result
                    if eq(combined.start,combined.end):
                        continue

                    # Existing one? => Remove it, add information with this combination
                    elif eq(res[-1][2].start,res[-1][2].end):
                        res.pop()
                    
                    # b) Multi item domains, remove ambiguity
                    else:
                        combined.left_open = not res[-1][2].right_open
            
            res.append(
                (
                    l.expr,
                    r.expr,
                    combined
                )
            )
    return res


