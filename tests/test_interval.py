#!/usr/bin/env python

"""Tests for `vodes.error.interval` class."""

from random import randrange
from vodes.symbolic.symbols import Boundary, BoundedExpression, BoundedValue, BoundedVariable, DummyVariable
from vodes.symbolic.power import Power
import pytest
from vodes.symbolic.interval import Interval
from pymbolic import var

# Utility functions for comparison
from tests.utils import assert_bounded_iv_equations

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


def test_interval_pow1(evaluators):
    #Arrange
    u1 = randrange(start=2,stop=22,step=2) #Even exponent
    x = var("x") 
    
    symbol = BoundedVariable(
        x.name,
        boundary=Boundary(
            lower=BoundedValue(value=-2,open=False),
            upper=BoundedValue(value=2,open=False)
        )
    )

    p = Interval(x-1,x+1) ** u1

    e = [
        BoundedExpression(
            boundary=Boundary(
                lower=BoundedValue(value=-2,open=False),
                upper=BoundedValue(value=0,open=True)
            ),
            expression=Interval(
                lower=Power(x+1,u1),
                upper=Power(x-1,u1)
            )
        ),
        BoundedExpression(
            boundary=Boundary(
                lower=BoundedValue(value=0,open=False),
                upper=BoundedValue(value=2,open=False)
            ),
            expression=Interval(
                lower=Power(x-1,u1),
                upper=Power(x+1,u1)
            )
        ),
    ]

    context = {}

    for eval in evaluators:
        #Act
        a = eval(context,symbol)(p)

        #Assert
        assert_bounded_iv_equations(a,e)

def test_interval_pow2(evaluators):
    #Arrange
    u1 = randrange(1,21,step=2) #Uneven exponent
    x = var("x") 
    
    symbol = BoundedVariable(
        x.name,
        boundary=Boundary(
            lower=BoundedValue(value=-2,open=False),
            upper=BoundedValue(value=2,open=False)
        )
    )

    p = Interval(x-1,x+1) ** u1

    e = [
        BoundedExpression(
            boundary=Boundary(
                lower=BoundedValue(value=-2,open=False),
                upper=BoundedValue(value=2,open=False)
            ),
            expression=Interval(
                lower=Power(x-1,u1),
                upper=Power(x+1,u1)
            )
        ),
    ]

    context = {}

    for eval in evaluators:
        #Act
        a = eval(context,symbol)(p)

        #Assert
        assert_bounded_iv_equations(a,e)
        
        
def test_interval_pow3(evaluators):
    #Arrange
    u1 = randrange(1,121)
    u2 = randrange(1,11,step=2) # Uneven exponent
    x = var("x") 
    
    symbol = BoundedVariable(
        x.name,
        boundary=Boundary(
            lower=BoundedValue(value=0,open=True),
            upper=BoundedValue(value=1,open=True)
        )
    )

    p = Interval(x)**u1 + Interval(-x,x)**u2

    e = [
        BoundedExpression(
            boundary=Boundary(
                lower=BoundedValue(value=0,open=True),
                upper=BoundedValue(value=1,open=True)
            ),
            expression=Interval(
                lower=x**u1 - x**u2,
                upper=x**u1 + x**u2
            )
        ),
    ]

    context = {}

    for eval in evaluators:
        #Act
        a = eval(context,symbol)(p)

        #Assert
        assert_bounded_iv_equations(a,e)


