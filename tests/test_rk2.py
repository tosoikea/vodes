#!/usr/bin/env python

"""Tests for `vodes` package."""

import pytest
import vodes.runge
import math

def test_rk2_f1():
    f = lambda t,y : 3 * y
    l = vodes.runge.RK2(f, 1, (0,3), {"dt":1})
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
    f = lambda t,y : 3 * y
    l = vodes.runge.RK2(f, 1, (0,3), {"dt":0.5})
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
    f = lambda t,y : t * y + t
    l = vodes.runge.RK2(f, 1, (0,5), {"dt":1})
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
    f = lambda t,y : t * y + t
    l = vodes.runge.RK2(f, 1, (0,3), {"dt":0.5})
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
