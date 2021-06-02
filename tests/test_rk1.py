#!/usr/bin/env python

"""Tests for `vodes` package."""

import pytest
from vodes.ode import runge
from sympy import symbols
from sympy.functions import exp

# define functions
y = symbols('y')
t = symbols('t')

f1 = 3 * y
f1_c = 3j * y
f2 = t * y + t

def test_rk1_f1():
    l = runge.RK1(f1, 1, (0,3), {"dt":1}, symbols=[y])
    l.calculate()

    actual = l.solution
    expected = [
        (0,1),
        (1,4),
        (2,16),
        (3,64)
    ]
    
    assert actual == expected

def test_rk1_f1_c():
    l = runge.RK1(f1_c, 1, (0,3), {"dt":1}, symbols=[y])
    l.calculate()

    actual = l.solution
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
    l = runge.RK1(f1, 1, (0,3), {"dt":0.5}, symbols=[y])
    l.calculate()

    actual = l.solution
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
    l = runge.RK1(f2, 1, (0,5), {"dt":1}, symbols=[y, t])
    l.calculate()

    actual = l.solution
    expected = [
        (0,1),
        (1,1),
        (2,3),
        (3,11),
        (4,47),
        (5,239)
    ]

    assert actual == expected

def test_rk1_v2():
    l = runge.RK1(f2, 1, (0,3), {"dt":0.5}, symbols=[y, t])
    l.calculate()

    actual = l.solution
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
