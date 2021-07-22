#!/usr/bin/env python

"""Tests for `vodes.error` package."""

from random import randrange
from typing import List
from vodes.symbolic.interval import Interval
from pymbolic.interop.sympy import PymbolicToSympyMapper

from pymbolic.primitives import Sum
from sympy import sqrt
from sympy.core.function import expand
from sympy.core.numbers import E
from sympy.core.symbol import symbols
import pytest
from vodes.symbolic.symbols import Boundary, BoundedValue, MachineError, BoundedExpression
from vodes.symbolic.interval import ExactIntersectionEvaluator as EIE, Interval
from pymbolic import var
from interval import interval, inf, imath


def __assert_bounds(actual, expected):
    assert len(actual) == len(expected)

    for i in range(len(actual)):
        assert isinstance(actual[i],BoundedExpression)
        assert actual[i].bound == expected[i]

def __assert_equations(actual,bounds,equations):
    __assert_bounds(actual, bounds)

    for i in range(len(actual)):
        assert expand(PymbolicToSympyMapper()(actual[i].expr.low)) == equations[i][0]
        assert expand(PymbolicToSympyMapper()(actual[i].expr.up)) == equations[i][1]


def test_multiple_intersections1():
    # Arrange
    e = MachineError()

    ##20x^3 - 20^x2 + 2
    p1 = 10*e**3 - 10*e**2 + 1
    ##20x^3 - 10x^2
    p2 = 10*e**3 - 5*e**2

    ## not a valid interval
    iv = Interval(p1,p2)
    expr = iv * 2

    # Act
    actual = EIE(context={}, symbol=MachineError())(expr)

    # Assert
    err = symbols('eps')
    intersection = sqrt(5)/5

    exp_equations = [
        [
            # [20e^3 - 10e^2, 20e^3 - 20e^2 + 2]
            20*err**3 - 10*err**2, 20*err**3 - 20*err**2 + 2
            ],
        [
            # [20e^3 - 20e^2 + 2, 20e^3 - 10e^2]
            20*err**3 - 20*err**2 + 2, 20*err**3-10*err**2
            ]
        ]

    exp_bounds = [
        Boundary(lower=MachineError().bound.lower,upper=BoundedValue(value=intersection,open=True)),
        Boundary(lower=BoundedValue(value=intersection, open=False),upper=MachineError().bound.upper),
        ]

    __assert_equations(actual=actual, bounds=exp_bounds, equations=exp_equations)