from vodes.symbolic.symbols import BoundedExpression
from pymbolic.interop.sympy import PymbolicToSympyMapper
from pymbolic.primitives import Expression
from sympy import expand


def __assert_equivalent(pl:Expression,pr:Expression):
    assert expand(PymbolicToSympyMapper()(pl)) == expand(PymbolicToSympyMapper()(pr))

def assert_bounds(actual, expected):
    assert len(actual) == len(expected)

    for i in range(len(actual)):
        assert isinstance(actual[i],BoundedExpression)
        assert actual[i].bound == expected[i]

def assert_bounded_iv_equations(actual,expected):
    assert len(actual) == len(expected)

    for i in range(len(actual)):
        # Equivalent Interval Expressions
        __assert_equivalent(actual[i].expr.low, expected[i].expr.low)
        __assert_equivalent(actual[i].expr.up, expected[i].expr.up)

        # Equivalent Bounds
        assert actual[i].bound == expected[i].bound