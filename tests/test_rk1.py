#!/usr/bin/env python

"""Tests for `vodes` package."""

import pytest
from vodes.ode.runge import RK1
from vodes.ode.problem import Problem
from pymbolic.primitives import Variable

# define functions
y = Variable('y')
t = Variable('t')

f1 = 3 * y
f1_c = 3j * y
f2 = t * y + t

def test_rk1_f1():
    p1 = Problem(f1, 1, (0,3))
    l = RK1(p1, {"dt":1})
    l.calculate()

    actual = l.solutions
    expected = [
        (0,1),
        (1,4),
        (2,16),
        (3,64)
    ]
    
    assert actual == expected

def test_rk1_f1_c():
    p1c = Problem(f1_c, 1, (0,3))
    l = RK1(p1c, {"dt":1})
    l.calculate()

    actual = l.solutions
    expected = [
        (0, 1), 
        (1, (1+3j)), 
        (2, (-8+6j)), 
        (3, (-26-18j))
    ]

    cactual = [(t,complex(y)) for (t,y) in actual]
    cexpected = [(t,complex(y)) for (t,y) in expected]
    
    assert cactual == cexpected

def test_rk1_v1():
    p1v1 = Problem(f1, 1, (0,3))
    l = RK1(p1v1, {"dt":0.5})
    l.calculate()

    actual = l.solutions
    expected = [
        (0,1),
        (0.5,2.5),
        (1,6.25),
        (1.5,15.625),
        (2,39.0625),
        (2.5,97.65625),
        (3,244.140625)
    ]
    
    assert actual == expected

def test_rk1_f2():
    p2 = Problem(f2, 1, (0,5))
    l = RK1(p2, {"dt":1})
    l.calculate()

    actual = l.solutions
    expected = [
        (0,1),
        (1,1),
        (2,3),
        (3,11),
        (4,47),
        (5,239)
    ]
    print(actual)
    assert actual == expected

def test_rk1_v2():
    p2v2 = Problem(f2, 1, (0,3))
    l = RK1(p2v2, {"dt":0.5})
    l.calculate()

    actual = l.solutions
    expected = [
        (0,1),
        (0.5,1),
        (1,1.5),
        (1.5,2.75),
        (2,5.5625),
        (2.5,12.125),
        (3,28.53125)
    ]
    
    assert actual == expected
