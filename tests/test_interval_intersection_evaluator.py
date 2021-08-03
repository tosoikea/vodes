#!/usr/bin/env python

"""Tests for `vodes.error` package."""
import pytest

from random import randrange
from typing import List
from vodes.symbolic.interval import Interval

from sympy import sqrt
from vodes.symbolic.symbols import Boundary, BoundedValue, MachineError, BoundedExpression
from vodes.symbolic.interval import Interval
from vodes.symbolic.mapper.exact_intersection_evaluator import ExactIntersectionEvaluator as EIE

# Utility functions for comparison
from tests.utils import assert_bounded_iv_equations


def test_multiple_intersections1():
    # Arrange
    e = MachineError()

    ##10x^3 - 10x^2 + 1
    p1 = 10*e**3 - 10*e**2 + 1
    ##10x^3 - 5x^2
    p2 = 10*e**3 - 5*e**2

    ## not a valid interval
    iv = Interval(p1,p2)
    expr = iv * 2

    # Act
    actual = EIE(context={}, symbol=MachineError())(expr)

    # Assert
    intersection = sqrt(5)/5

    expected = [
        BoundedExpression(
            # [20e^3 - 10e^2, 20e^3 - 20e^2 + 2]
            expression=Interval(
                lower=20*e**3 - 10*e**2,
                upper=20*e**3 - 20*e**2 + 2
            ),
            boundary=Boundary(
                lower=MachineError().bound.lower,
                upper=BoundedValue(value=intersection,open=True)
            )
        ),
        BoundedExpression(
            # [20e^3 - 20e^2 + 2, 20e^3 - 10e^2]
            expression=Interval(
               lower=20*e**3 - 20*e**2 + 2,
               upper= 20*e**3-10*e**2
            ),
            boundary=Boundary(
                lower=BoundedValue(value=intersection, open=False),
                upper=MachineError().bound.upper
            )
        )
    ]

    assert_bounded_iv_equations(actual=actual, expected=expected)