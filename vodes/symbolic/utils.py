from pymbolic.primitives import Expression

def gt(o,v) -> bool:
    if isinstance(o,Expression):
        from pymbolic import evaluate
        return evaluate(o.gt(v))
    elif isinstance(v,Expression):
        from pymbolic import evaluate
        return evaluate(v.lt(o))
    else:
        return o > v

def ge(o,v) -> bool:
    if isinstance(o,Expression):
        from pymbolic import evaluate
        return evaluate(o.ge(v))
    elif isinstance(v,Expression):
        from pymbolic import evaluate
        return evaluate(v.le(o))
    else:
        return o >= v

def lt(o,v) -> bool:
    if isinstance(o,Expression):
        from pymbolic import evaluate
        return evaluate(o.lt(v))
    elif isinstance(v,Expression):
        from pymbolic import evaluate
        return evaluate(v.gt(o))
    else:
        return o < v

def le(o,v) -> bool:
    if isinstance(o,Expression):
        from pymbolic import evaluate
        return evaluate(o.le(v))
    elif isinstance(v,Expression):
        from pymbolic import evaluate
        return evaluate(v.ge(o))
    else:
        return o <= v

def eq(o,v) -> bool:
    if isinstance(o,Expression):
        from pymbolic import evaluate
        return evaluate(o.eq(v))
    elif isinstance(v,Expression):
        from pymbolic import evaluate
        return evaluate(v.eq(o))
    else:
        return o == v