#!/usr/bin/env python

"""Tests for `vodes` package."""

import pytest
import vodes.runge

def test_rk1_f1():
    f = lambda t,y : 3 * y
    l = vodes.runge.RK1(f, 1, (0,3), {"dt":1})
    l.calculate()

    actual = l.solution
    expected = [
        (0,1),
        (1,4),
        (2,16),
        (3,64)
    ]
    
    assert actual == expected

def test_rk1_v1():
    f = lambda t,y : 3 * y
    l = vodes.runge.RK1(f, 1, (0,3), {"dt":0.5})
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
    f = lambda t,y : t * y + t
    l = vodes.runge.RK1(f, 1, (0,5), {"dt":1})
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
    f = lambda t,y : t * y + t
    l = vodes.runge.RK1(f, 1, (0,3), {"dt":0.5})
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
