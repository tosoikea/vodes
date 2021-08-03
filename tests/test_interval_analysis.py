#!/usr/bin/env python
import pytest
"""Tests for `vodes.error` package."""

import math
from random import randrange

from sympy.core.power import Pow
from vodes.symbolic.interval import Interval

from pymbolic.primitives import Power, Quotient
from sympy.core.symbol import symbols
from vodes.symbolic.symbols import Boundary, BoundedExpression, BoundedValue, MachineError
from vodes.error.analysis import IntervalAnalysis
from pymbolic import var
from interval import interval

from pymbolic.primitives import Expression


def __assert_expressions(l,r):
    from pymbolic.interop.sympy import PymbolicToSympyMapper
    from sympy import simplify

    if isinstance(l,Expression):
        l = PymbolicToSympyMapper()(l)
    
    if isinstance(r, Expression):
        r = PymbolicToSympyMapper()(r)

    assert simplify(l-r) == 0

def __bound(iv):
    b = max(abs(iv))
    return b[0]

def __assert_constant(actual,expected):
    assert(len(actual)==1)
    assert(isinstance(actual[0],BoundedExpression))

    math.isclose(actual[0].expr, __bound(expected), rel_tol=0.01)

def __assert_equation(actual,expected):
    assert len(actual) == 1
    assert isinstance(actual[0],BoundedExpression)

    __assert_expressions(actual[0].expr, expected)

def __assert_bounds(actual, expected):
    assert len(actual) == len(expected)

    for i in range(len(actual)):
        assert isinstance(actual[i],BoundedExpression)
        __assert_expressions(actual[i].bound.lower.value, expected[i].lower.value)
        __assert_expressions(actual[i].bound.upper.value, expected[i].upper.value)

        assert actual[i].bound.lower.open == expected[i].lower.open
        assert actual[i].bound.upper.open == expected[i].upper.open

def __assert_equations(actual,bounds,equations):
    __assert_bounds(actual, bounds)

    for i in range(len(actual)):
        __assert_expressions(actual[i].expr, equations[i])

def test_iv_full_sub_1():
    """Comparison of interval analysis with direct usage of interval library"""

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
    __assert_constant(actual,expected)


def test_iv_symbolic_1():
    # Arrange
    x = var("x")
    p = x * 2
    xv = 5

    # Act
    ia = IntervalAnalysis(p)
    actual = ia.absolute(context={
        "x": xv
    },
    min_precision=1,
    max_precision=113)

    err = symbols('eps')
    # 20e + 10e^2
    expected = 20*err + 10*err**2

    __assert_equation(actual, expected)
    
def test_iv_symbolic_2():
    # Arrange
    x = var("x")
    p = x + 20
    xv = 10

    # Act
    ia = IntervalAnalysis(p)
    actual = ia.absolute(context={
        "x": xv
    },
    min_precision=1,
    max_precision=113)

    err = symbols('eps')
    # 40e + 10e^2
    expected = 40*err + 10*err**2

    __assert_equation(actual, expected)


def test_iv_symbolic_3():
    # Arrange
    x = var("x")
    p = x**2
    xv = 5

    # Act
    ia = IntervalAnalysis(p)
    actual = ia.absolute(context={
        "x": xv
    },
    min_precision=1,
    max_precision=113)

    err = symbols('eps')
    # 25e^3+75e^2+75e
    expected = 25*err**3 + 75*err**2 + 75*err

    __assert_equation(actual, expected)

def test_iv_symbolic_4():
    # Arrange
    x = var("x")
    p = (x-1)
    xv = Quotient(3,2)

    # Act
    ia = IntervalAnalysis(p)
    actual = ia.absolute(context={
        "x": xv
    },
    min_precision=1,
    max_precision=113)

    # Assert

    #
    b1 = Quotient(1,3)
    exp_bounds = [
        Boundary(lower=MachineError(max_precision=113).bound.lower,upper=BoundedValue(value=b1,open=True)),
        Boundary(lower=BoundedValue(value=b1, open=False),upper=MachineError(min_precision=1).bound.upper),
        ]

    err = symbols('eps')
    exp_equations = [
        # - 1/2 + (1.5e + 0.5)(1+e)
        (-1) * Quotient(1,2) + (Quotient(3,2) * err + Quotient(1,2)) * (1 + err),
        (-1) * Quotient(1,2) + (Quotient(3,2) * err + Quotient(1,2)) * (1 + err)
    ]

    __assert_equations(actual, bounds=exp_bounds, equations=exp_equations)

def test_iv_symbolic_5():
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
    b1 = Pow(3,-1)

    exp_bounds = [
        Boundary(
            lower=MachineError(max_precision=max_prec).bound.lower,
            upper=BoundedValue(value=b1,open=True)
        ),
        Boundary(
            lower=BoundedValue(value=b1,open=False),
            upper=MachineError(min_precision=min_prec).bound.upper
        )
    ]

    err = symbols('eps')
    exp_equations = [
        # -0.25 + (1.5e**2+2*e+0.5)**2*(1+e)
        (-1) * Quotient(1,4) + (Quotient(3,2)*err**2 + 2*err + Quotient(1,2))**2*(1+err),
        (-1) * Quotient(1,4) + (Quotient(3,2)*err**2 + 2*err + Quotient(1,2))**2*(1+err)
    ]

    __assert_equations(actual, bounds=exp_bounds, equations=exp_equations)


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
            Boundary(
                lower=BoundedValue(value=Power(2,-max_prec),open=False),
                upper=BoundedValue(value=Power(2,-min_prec),open=False)
            )
            ]
        )

def test_iv_symbolic_7():
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
    exp_bounds = [
        Boundary(lower=MachineError(max_precision=max_prec).bound.lower,upper=MachineError(min_precision=min_prec).bound.upper),
        ]
    exp_equations = [0]

    __assert_equations(actual, bounds=exp_bounds, equations=exp_equations)