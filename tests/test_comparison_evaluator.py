#!/usr/bin/env python

"""Tests for `vodes.error.interval` class."""
from math import exp
import pytest

from random import randrange

from vodes.symbolic.mapper.comparison_evaluator import ComparisonEvaluator as CE

# Custom Expression Library
from vodes.symbolic.expressions.bounded import Domain, BoundedExpression, BoundedVariable, DummyVariable, MachineError
from vodes.symbolic.expressions.interval import Interval
from vodes.symbolic.expressions.nthroot import NthRoot

# Symbolic Expression Library
from pymbolic import var
from pymbolic.primitives import Quotient, Power

# Utility functions for comparison
from tests.utils import assert_bounded_iv_equations

def test_minimum1():
    # minimum(x,x,x,x)

    # Arrange
    d = Domain(-1,1)
    x = BoundedVariable('x',d)
    exprs = [
        x,x,x,x
    ]

    # Act
    res = CE({},x)._minimum(exprs=exprs,boundary=d)

    # Assert
    expected = [
        BoundedExpression(exprs[0],d)
    ]

    assert expected == res

def test_minimum2():
    # minimum( (x-1)^2, (x+1)^2 )
    from vodes.symbolic.utils import sort_bounded

    # Arrange
    d = Domain(-2,2)
    x = BoundedVariable('x',d)
    exprs = [
        (x-1)**2,(x+1)**2
    ]

    # Act
    res = CE({},x)._minimum(exprs=exprs,boundary=d)

    # Assert
    expected = [
        BoundedExpression(exprs[1],Domain(-2,0,False,True)),
        BoundedExpression(exprs[0],Domain(0,2,False,False))
    ]

    sort_bounded(res)
    assert expected == res

def test_contains1():
    # Arrange
    d = Domain(-2,2)
    x = BoundedVariable('x',d)

    expr = Interval(0,2)
    val = 1

    # Act
    res = CE({},x).contains(iv=expr,xs=val,d=d)

    # Assert
    expected = [
        (d,[1])
    ]

    assert expected == res
    
def test_contains2():
    # Arrange
    d = Domain(-2,2)
    x = BoundedVariable('x',d)

    expr = Interval(x-1,x+1)
    val = 0

    # Act
    res = CE({},x).contains(iv=expr,xs=val,d=d)

    # Assert
    expected = [
        (Domain(-2,-1,right_open=True),[]),
        (Domain(-1,1),[0]),
        (Domain(1,2,left_open=True),[])
    ]

    assert expected == res

def test_contains3():
    from vodes.symbolic.expressions.constants import Pi

    # Arrange
    d = Domain(-200,200)
    x = BoundedVariable('x',d)

    expr = Interval(0,x)
    val = Pi()

    # Act
    res = CE({},x).contains(iv=expr,xs=val,d=d)

    # Assert
    expected = [
        (Domain(-200,Pi(),right_open=True),[]),
        (Domain(Pi(),200),[Pi()])
    ]

    assert expected == res

def test_multiple_intersections1():
    # Arrange
    min_prec = 1
    max_prec = 113
    e = MachineError(min_precision=min_prec,max_precision=max_prec)

    ##10x^3 - 10x^2 + 1
    p1 = 10*e**3 - 10*e**2 + 1
    ##10x^3 - 5x^2
    p2 = 10*e**3 - 5*e**2

    ## not a valid interval
    iv = Interval(p1,p2)
    expr = iv * 2

    # Act
    actual = CE(context={}, symbol=e)(expr)

    # Assert
    intersection = Quotient(Power(5,Quotient(1,2)),5)

    expected = [
        BoundedExpression(
            # [20e^3 - 10e^2, 20e^3 - 20e^2 + 2]
            expression=Interval(
                lower=20*e**3 - 10*e**2,
                upper=20*e**3 - 20*e**2 + 2
            ),
            boundary=Domain(
                start=e.bound.start,
                left_open=e.bound.left_open,
                end=intersection
            )
        ),
        BoundedExpression(
            # [20e^3 - 20e^2 + 2, 20e^3 - 10e^2]
            expression=Interval(
               lower=20*e**3 - 20*e**2 + 2,
               upper= 20*e**3-10*e**2
            ),
            boundary=Domain(
                start=intersection,
                end=e.bound.end,
                left_open=True,
                right_open=e.bound.right_open
            )
        )
    ]

    print(list(map(str,actual)))
    print(len(actual))
    print(len(expected))

    assert_bounded_iv_equations(actual=actual, expected=expected)