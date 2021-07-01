#!/usr/bin/env python

"""Tests for `vodes.error.interval` class."""

from random import randrange
from vodes.symbolic.interval import Interval, IntervalEvaluator as IE


def test_setup():
    # Arrange
    upper = randrange(5000)
    lower = upper - randrange(5000)

    # Act 
    i = Interval(lower,upper)

    # Assert
    assert i.upper == upper
    assert i.lower == lower

def test_scalar_add():
    # Arrange
    upper = randrange(5000)
    lower = upper - randrange(5000)
    a = randrange(5000)

    # Act
    li = IE()(Interval(lower,upper) + a)
    ri = IE()(a + Interval(lower,upper))

    # Assert
    assert li.upper == upper + a
    assert ri.upper == upper + a
    assert li.lower == lower + a
    assert ri.lower == lower + a

def test_scalar_sub():
    # Arrange
    upper = randrange(5000)
    lower = upper - randrange(5000)
    a = randrange(5000)

    # Act
    li = IE()(Interval(lower,upper) - a)
    ri = IE()(a - Interval(lower,upper))

    # Assert
    assert li.upper == upper - a
    assert ri.upper == a - lower
    assert li.lower == lower - a
    assert ri.lower == a - upper

def test_scalar_mul():
    # Arrange
    upper = randrange(5000)
    lower = upper - randrange(5000)
    a = randrange(5000)

    # Act
    li = IE()(Interval(lower,upper) * a)
    ri = IE()(a * Interval(lower,upper))

    # Assert
    assert li.upper == upper * a
    assert ri.upper == upper * a
    assert li.lower == lower * a
    assert ri.lower == lower * a

def test_interval_add():
    # Arrange
    u1 = randrange(5000)
    l1 = u1 - randrange(5000)

    u2 = randrange(5000)
    l2 = u2 - randrange(5000)

    # Act
    li = IE()(Interval(l1,u1) + Interval(l2,u2))
    ri = IE()(Interval(l2,u2) + Interval(l1,u1))

    # Assert
    assert li.upper == u1 + u2
    assert ri.upper == u1 + u2
    assert li.lower == l1 + l2
    assert ri.lower == l1 + l2
    
def test_interval_sub():
    # Arrange
    u1 = randrange(5000)
    l1 = u1 - randrange(5000)

    u2 = randrange(5000)
    l2 = u2 - randrange(5000)

    # Act
    li = IE()(Interval(l1,u1) - Interval(l2,u2))
    ri = IE()(Interval(l2,u2) - Interval(l1,u1))

    # Assert
    assert li.upper == u1 - l2
    assert ri.upper == u2 - l1
    assert li.lower == l1 - u2
    assert ri.lower == l2 - u1
    
    
def test_interval_mul():
    # Arrange
    u1 = randrange(5000)
    l1 = u1 - randrange(5000)

    u2 = randrange(5000)
    l2 = u2 - randrange(5000)

    # Act
    li = IE()(Interval(l1,u1) * Interval(l2,u2))
    ri = IE()(Interval(l2,u2) * Interval(l1,u1))

    # Assert
    assert li.upper == max(u1 * l2, u1 * u2, l1 * l2, l1 * u2)
    assert ri.upper == max(u1 * l2, u1 * u2, l1 * l2, l1 * u2)
    assert li.lower == min(u1 * l2, u1 * u2, l1 * l2, l1 * u2)
    assert ri.lower == min(u1 * l2, u1 * u2, l1 * l2, l1 * u2)

