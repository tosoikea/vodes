#!/usr/bin/env python

"""Tests for `vodes.runge` package."""

import pytest
from vodes.ode.runge import RK2
from vodes.ode.problem import Problem
import math
from sympy import symbols
from sympy.functions import exp

# define functions
y = symbols('y')
t = symbols('t')

f1 = 3 * y
f1_c = 3j * y
f2 = t * y + t


def test_rk2_f1():
    p1 = Problem(f1, 1, (0,3))
    l = RK2(p1, {"dt":1})
    l.calculate()

    actual = l.solution
    expected = [
        (0,1),
        (1,8.5),
        (2,72.25),
        (3,614.125)
    ]
    
    assert actual == expected

def test_rk2_v1():
    p1 = Problem(f1, 1, (0,3))
    l = RK2(p1, {"dt":0.5})
    l.calculate()

    actual = l.solution
    expected = [
        (0,1),
        (0.5,3.625),
        (1,13.140625),
        (1.5,47.634765625),
        (2,172.676025390625),
        (2.5,625.950592041016),
        (3,2269.07089614868)
    ]
    
    assert len(actual) == len(expected)
    assert all([a[0] == b[0] and math.isclose(a[1],b[1],rel_tol=0.05) for a, b in zip(actual, expected)])

def test_rk2_f2():
    p2 = Problem(f2, 1, (0,5))
    l = RK2(p2, {"dt":1})
    l.calculate()

    actual = l.solution
    expected = [
        (0,1),
        (1,2),
        (2,9.5),
        (3,67.25),
        (4,715.625),
        (5,11106.6875)
    ]

    assert actual == expected

def test_rk2_v2():
    p2 = Problem(f2, 1, (0,3))
    l = RK2(p2, {"dt":0.5})
    l.calculate()

    actual = l.solution
    expected = [
        (0,1),
        (0.5,1.25),
        (1,2.234375),
        (1.5,4.8623046875),
        (2,12.190185546875),
        (2.5,35.2730102539063),
        (3,119.154346466064)
    ]
    
    assert len(actual) == len(expected)
    assert all([a[0] == b[0] and math.isclose(a[1],b[1],rel_tol=0.05) for a, b in zip(actual, expected)])
