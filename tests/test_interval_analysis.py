#!/usr/bin/env python
import pytest
"""Tests for `vodes.error` package."""

import math
from random import randrange

from vodes.error.analysis import IntervalAnalysis
from interval import interval

# Custom Expression Library
from vodes.symbolic.expressions.bounded import Domain, BoundedExpression, MachineError

# Symbolic Expression Library
from pymbolic import var
from pymbolic.primitives import Power, Quotient

from sympy.core.power import Pow

def __assert_bounds(actual, expected):
    from tests.utils import assert_bounded_equations

    actual_noexpr = []
    for i in range(len(actual)):
        assert isinstance(actual[i],BoundedExpression)
        actual_noexpr.append(
            BoundedExpression(
                expression=0,
                boundary=actual[i].bound
            )
        )
    
    expected_noexpr = [BoundedExpression(expression=0,boundary=bound) for bound in expected]
    assert_bounded_equations(actual_noexpr,expected_noexpr)

def test_iv_full_sub_1():
    """Comparison of interval analysis with direct usage of interval library"""
    from tests.utils import assert_constant

    # Arrange
    x = var("x")
    p = x * randrange(20)
    xv = randrange(5)
    ev = 2**-(randrange(1,10))

    # Act
    ia = IntervalAnalysis(p)
    actual = ia.absolute(context={
        "x": xv,
        MachineError().name: ev
    },
    min_precision=1,
    max_precision=113)
    
    expr_p_iv = (interval(xv,xv) * 2) 
    expr_err_iv = (interval(xv,xv) * interval(1-ev,1+ev)) * 2 * interval(1-ev,1+ev)
    expected = expr_p_iv - expr_err_iv

    # Assert
    assert_constant(actual,expected)


def test_iv_symbolic_1():
    from tests.utils import assert_bounded_equations

    # Arrange
    x = var("x")
    p = x * 2
    xv = 5

    min_prec = 1
    max_prec = 113

    # Act
    ia = IntervalAnalysis(p)
    actual = ia.absolute(context={
        "x": xv
    },
    min_precision=min_prec,
    max_precision=max_prec)

    err = MachineError(min_precision=min_prec,max_precision=max_prec)

    # 20e + 10e^2
    expected = [
        BoundedExpression(
            expression=20*err + 10*err**2,
            boundary=Domain(
                start=err.bound.start,
                left_open=err.bound.left_open,
                end=err.bound.end,
                right_open=err.bound.right_open
            )
        )
    ]

    assert_bounded_equations(actual, expected)
    
def test_iv_symbolic_2():
    from tests.utils import assert_bounded_equations

    # Arrange
    x = var("x")
    p = x + 20
    xv = 10

    min_prec = 1
    max_prec = 113

    # Act
    ia = IntervalAnalysis(p)
    actual = ia.absolute(context={
        "x": xv
    },
    min_precision=min_prec,
    max_precision=max_prec)

    err = MachineError(min_precision=min_prec,max_precision=max_prec)

    # 40e + 10e^2
    expected = [
        BoundedExpression(
            expression=40*err + 10*err**2,
            boundary=Domain(
                start=err.bound.start,
                left_open=err.bound.left_open,
                end=err.bound.end,
                right_open=err.bound.right_open
            )
        )
    ]

    assert_bounded_equations(actual, expected)


def test_iv_symbolic_3():
    from tests.utils import assert_bounded_equations

    # Arrange
    x = var("x")
    p = x**2
    xv = 5

    min_prec = 1
    max_prec = 113

    # Act
    ia = IntervalAnalysis(p)
    actual = ia.absolute(context={
        "x": xv
    },
    min_precision=min_prec,
    max_precision=max_prec)

    err = MachineError(min_precision=min_prec,max_precision=max_prec)

    # 25e^3+75e^2+75e
    expected = [
        BoundedExpression(
            expression=25*err**3 + 75*err**2 + 75*err,
            boundary=Domain(
                start=err.bound.start,
                left_open=err.bound.left_open,
                end=err.bound.end,
                right_open=err.bound.right_open
            )
        )
    ]

    assert_bounded_equations(actual, expected)

def test_iv_symbolic_4():
    from tests.utils import assert_bounded_equations

    # Arrange
    x = var("x")
    p = (x-1)
    xv = Quotient(3,2)

    min_prec = 1
    max_prec = 113

    # Act
    ia = IntervalAnalysis(p)
    actual = ia.absolute(context={
        "x": xv
    },
    min_precision=min_prec,
    max_precision=max_prec)

    # Assert
    err = MachineError(min_precision=min_prec,max_precision=max_prec)
    b1 = Quotient(1,3)

    # - 1/2 + (1.5e + 0.5)(1+e)
    expected = [
        BoundedExpression(
            expression=(-1) * Quotient(1,2) + (Quotient(3,2) * err + Quotient(1,2)) * (1 + err),
            boundary=Domain(
                start=err.bound.start,
                left_open=err.bound.left_open,
                end=b1,
                right_open=True
            ),
        ),
        BoundedExpression(
            expression=(-1) * Quotient(1,2) + (Quotient(3,2) * err + Quotient(1,2)) * (1 + err),
            boundary=Domain(
                start=b1,
                end=err.bound.end,
                right_open=err.bound.right_open
            ),
        )
    ]

    assert_bounded_equations(actual, expected)

def test_iv_symbolic_5():
    from tests.utils import assert_bounded_equations
    
    # Arrange
    x = var("x")
    p = (x-1)**2
    xv = Quotient(3,2)

    max_prec = 113
    min_prec = 1

    # Act
    ia = IntervalAnalysis(p)
    actual = ia.absolute(context={
        "x": xv
    },
    min_precision=min_prec,
    max_precision=max_prec)

    # Assert
    err = MachineError(min_precision=min_prec,max_precision=max_prec)
    b1 = Pow(3,-1)

    expected = [
        BoundedExpression(
            expression=(-1) * Quotient(1,4) + (Quotient(3,2)*err**2 + 2*err + Quotient(1,2))**2*(1+err),
            boundary=Domain(
                start=err.bound.start,
                left_open=err.bound.left_open,
                end=b1,
                right_open=True
            )
        ),
        BoundedExpression(
            expression=(-1) * Quotient(1,4) + (Quotient(3,2)*err**2 + 2*err + Quotient(1,2))**2*(1+err),
            boundary=Domain(
                start=b1,
                end=err.bound.end,
                right_open=err.bound.right_open
            )
        )
    ]

    assert_bounded_equations(actual, expected)


def test_iv_symbolic_6():
    """Test case to validate, that the correct boundaries are returned for problems"""
    # Arrange
    x = var("x")
    p = x + 1

    max_prec = 113
    min_prec = 11

    xv = randrange(20)
    
    # Act
    ia = IntervalAnalysis(p)

    actual = ia.absolute(
        context={
            "x": xv
        },
        min_precision=min_prec,
        max_precision=max_prec
        )

    # Assert
    __assert_bounds(
        actual=actual,
        expected=[
            Domain(
                start=Power(2,-max_prec),
                end=Power(2,-min_prec)
            )
        ]
    )

def test_iv_symbolic_7():
    from tests.utils import assert_bounded_equations
    
    """Test case to validate edge case, where expression is constant"""
    # Arrange
    p = 0

    max_prec = 113
    min_prec = 11

    # Act
    ia = IntervalAnalysis(p)

    actual = ia.absolute(
        context={},
        min_precision=min_prec,
        max_precision=max_prec
        )    
        
    # Assert
    expected = [
        BoundedExpression(
            expression=0,
            boundary=Domain(
                start=MachineError(max_precision=max_prec).bound.start,
                left_open=MachineError(max_precision=max_prec).bound.left_open,
                end=MachineError(min_precision=min_prec).bound.end,
                right_open=MachineError(min_precision=min_prec).bound.right_open,
            )
        )
    ]

    assert_bounded_equations(actual, expected)