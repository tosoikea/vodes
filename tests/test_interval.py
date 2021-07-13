#!/usr/bin/env python

"""Tests for `vodes.error.interval` class."""

from random import randrange
from vodes.symbolic.symbols import BoundedExpression, BoundedVariable, DummyVariable
import pytest
from vodes.symbolic.interval import Interval

@pytest.fixture
def evaluators():
    from vodes.symbolic.interval import ExactIntersectionEvaluator as EIE
    return [
        lambda c,s : EIE(context=c, symbol=s)
    ]

@pytest.fixture
def static_evaluators():
    from vodes.symbolic.interval import ExactIntersectionEvaluator as EIE
    return [
       EIE(context={}, symbol=DummyVariable())
    ]

def assert_static(res,val):
    assert(len(res)==1)
    assert(isinstance(res[0],BoundedExpression))

    assert res[0].expr == val

def test_setup():
    # Arrange
    upper = randrange(5000)
    lower = upper - randrange(5000)

    # Act 
    i = Interval(lower,upper)

    # Assert
    assert i.up == upper
    assert i.low == lower

def test_scalar_add(static_evaluators):
    # Arrange
    upper = randrange(5000)
    lower = upper - randrange(5000)
    a = randrange(5000)

    e = Interval(lower=lower + a, upper=upper + a)

    for eval in static_evaluators:
        # Act
        li = eval(Interval(lower,upper) + a)
        ri = eval(a + Interval(lower,upper))

        # Assert
        assert_static(li, e)
        assert_static(ri, e)

def test_scalar_sub(static_evaluators):
    # Arrange
    upper = randrange(5000)
    lower = upper - randrange(5000)
    a = randrange(5000)

    el = Interval(lower = lower - a, upper = upper - a)
    er = Interval(lower = a - upper, upper = a - lower)

    for eval in static_evaluators:
        # Act
        li = eval(Interval(lower,upper) - a)
        ri = eval(a - Interval(lower,upper))

        # Assert
        assert_static(li, el)
        assert_static(ri, er)

def test_scalar_mul(static_evaluators):
    # Arrange
    upper = randrange(5000)
    lower = upper - randrange(5000)
    a = randrange(5000)

    e = Interval(lower = lower * a, upper = upper * a)

    for eval in static_evaluators:
        # Act
        li = eval(Interval(lower,upper) * a)
        ri = eval(a * Interval(lower,upper))

        # Assert
        assert_static(li, e)
        assert_static(ri, e)

def test_interval_add(static_evaluators):
    # Arrange
    u1 = randrange(5000)
    l1 = u1 - randrange(5000)

    u2 = randrange(5000)
    l2 = u2 - randrange(5000)

    for eval in static_evaluators:
        # Act
        li = eval(Interval(l1,u1) + Interval(l2,u2))
        ri = eval(Interval(l2,u2) + Interval(l1,u1))

        # Assert
        assert_static(li, Interval(lower=l1+l2,upper=u1+u2))
        assert_static(ri, Interval(lower=l1+l2,upper=u1+u2))
    
def test_interval_sub(static_evaluators):
    # Arrange
    u1 = randrange(5000)
    l1 = u1 - randrange(5000)

    u2 = randrange(5000)
    l2 = u2 - randrange(5000)

    for eval in static_evaluators:
        # Act
        li = eval(Interval(l1,u1) - Interval(l2,u2))
        ri = eval(Interval(l2,u2) - Interval(l1,u1))

        # Assert
        assert_static(li, Interval(lower=l1-u2,upper=u1-l2))
        assert_static(ri, Interval(lower=l2-u1,upper=u2-l1))
    
def test_interval_mul(static_evaluators):
    # Arrange
    u1 = randrange(5000)
    l1 = u1 - randrange(5000)

    u2 = randrange(5000)
    l2 = u2 - randrange(5000)

    e = Interval(
        lower=min(u1 * l2, u1 * u2, l1 * l2, l1 * u2),
        upper=max(u1 * l2, u1 * u2, l1 * l2, l1 * u2)
    )

    for eval in static_evaluators:
        # Act
        li = eval(Interval(l1,u1) * Interval(l2,u2))
        ri = eval(Interval(l2,u2) * Interval(l1,u1))

        # Assert
        assert_static(li, e)
        assert_static(ri, e)



