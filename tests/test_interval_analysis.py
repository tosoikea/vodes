#!/usr/bin/env python

"""Tests for `vodes.error` package."""

import math
from random import randrange
from vodes.symbolic.interval import Interval
from pymbolic.interop.sympy import PymbolicToSympyMapper

from pymbolic.primitives import Sum
from sympy.core.function import expand
from sympy.core.numbers import E
from sympy.core.symbol import symbols
from vodes.symbolic.symbols import BoundedExpression
import pytest
from vodes.error.analysis import IntervalAnalysis
from pymbolic import var
from interval import interval, inf, imath

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

    actual_sym = PymbolicToSympyMapper()(actual[0].expr)
    
    assert expand(actual_sym) == expected

def test_iv_full_sub_1():
    # Arrange
    x = var("x")
    p = x * randrange(20)
    xv = randrange(5)
    ev = 2**-(randrange(1,10))

    # Act
    ia = IntervalAnalysis(p)
    actual = ia.absolute(context={
        "x": xv,
        "e": ev
    })
    
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
    })

    err = symbols('e')
    # 20e + 10e²
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
    })

    err = symbols('e')
    # 40e + 10e²
    expected = 40*err + 10*err**2

    __assert_equation(actual, expected)

