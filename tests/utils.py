from vodes.symbolic.expressions.interval import Interval
from vodes.symbolic.expressions.bounded import BoundedExpression, Domain
from pymbolic.interop.sympy import PymbolicToSympyMapper
from pymbolic.primitives import Expression
from sympy import expand


def __assert_equivalent(pl:Expression,pr:Expression):
    if isinstance(pl,Expression):
        pl = expand(PymbolicToSympyMapper()(pl))

    if isinstance(pr,Expression):
        pr = expand(PymbolicToSympyMapper()(pr))

    assert pl == pr

def assert_bounds_equal(l:Domain,r:Domain):
    assert l.left_open == r.left_open
    assert l.right_open == r.right_open

    __assert_equivalent(l.start, r.start)
    __assert_equivalent(l.end, r.end)


def assert_bounded_iv_equations(actual,expected):
    assert len(actual) == len(expected)

    for i in range(len(actual)):
        assert isinstance(actual[i],BoundedExpression)
        assert isinstance(actual[i].expr,Interval)

        # Equivalent Interval Expressions
        __assert_equivalent(actual[i].expr.low, expected[i].expr.low)
        __assert_equivalent(actual[i].expr.up, expected[i].expr.up)

        # Equivalent Bounds
        assert_bounds_equal(actual[i].bound, expected[i].bound)

def assert_bounded_equations(actual,expected):
    assert len(actual) == len(expected)

    for i in range(len(actual)):
        assert isinstance(actual[i],BoundedExpression)

        # Equivalent Expressions
        __assert_equivalent(actual[i].expr, expected[i].expr)

        # Equivalent Bounds
        assert_bounds_equal(actual[i].bound, expected[i].bound)

def bound(iv):
    b = max(abs(iv))
    return b[0]

def assert_constant(actual,expected,rel_tol=0.01):
    import math

    assert(len(actual)==1)
    assert(isinstance(actual[0],BoundedExpression))

    math.isclose(actual[0].expr, bound(expected), rel_tol=rel_tol)