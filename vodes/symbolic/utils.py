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

def _sort_bounded_tuple(ls:list):
    """Sort a list of bounded expressions"""
    from functools import cmp_to_key
    list.sort(ls,key=cmp_to_key(lambda item1,item2: compare(item1[1],item2[1])))

def merge_unary(ls:list) -> List[Tuple[Expression,object]]:
    """Merges a single list of bounded expressions to tuples made up of the expression and the corresponding domain."""
    return _merge_unary([(l.expr,l.bound) for l in ls])

def _merge_unary(ls:list) -> List[Tuple[Expression,object]]:
    """Merges a single list of expressions associated with boundaries to tuples made up of the expression and the corresponding domain."""
    res = []

    _sort_bounded_tuple(ls)

    for (lexpr,lbound) in ls:
        if len(res) == 0:
            res.append(
                (lexpr,lbound)
            )
        else:
            (rexpr,rbound) = res.pop()
            existing = rbound.intersect(lbound)

            if not(existing is None):
                # expand existing boundary
                expansion = rbound.union(existing)
                # reduce new boundary
                reduction = lbound.difference(expansion)

                if len(expansion) != 1:
                    raise ValueError(f"Invalid expansion {expansion}")

                res.append(
                    (
                        rexpr,
                        expansion[0]
                    )
                )

                # No more boundary left
                if len(reduction) == 0:
                    continue
                elif len(reduction) == 1:
                    lbound = reduction[0]
                else:
                    raise ValueError(f"Invalid combination {reduction} (from exp. {expansion})")
            else:
                res.append((rexpr,rbound))

            res.append(
                (
                    lexpr,
                    lbound
                )
            )

    return res

def merge(ls:list,rs:list) -> List[Tuple[Tuple[Expression,Expression],object]]:
    """Merges two list of bounded expressions to tuples made up of a pair of items from the lists and the shared domain."""

    res = []
    mls = merge_unary(ls)
    mrs = merge_unary(rs)
    
    for (lexpr,lbound) in mls:
        for (rexpr,rbound) in mrs:
            combined = lbound.intersect(rbound)
            if combined is None:
                continue

            res.append(
                ((lexpr,rexpr),combined)
            )

    return _merge_unary(res)
