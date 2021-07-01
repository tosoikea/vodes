#!/usr/bin/env python

"""Tests for `vodes.error` package."""

import pytest
from vodes.error.analysis import IntervalAnalysis
from pymbolic import var
from interval import interval, inf, imath

def __bound(iv):
    b = max(abs(iv))
    return b[0]

def test_iv_full_sub_1():
    # Arrange
    x = var("x")
    p = x * 2
    xv = 5
    ev = 0.25

    # Act
    ia = IntervalAnalysis(p)
    actual = ia.absolute(context={
        "x": xv,
        "e": ev
    })
    
    expr_p_iv = (interval(xv,xv) * 2) 
    expr_err_iv = (interval(xv,xv) * interval(1-ev,1+ev)) * 2 * interval(1-ev,1+ev)
    expected = expr_p_iv - expr_err_iv

    # Assert
    assert actual == __bound(expected)


def test_iv_symbolic_1():
    # Arrange
    x = var("x")
    p = x * 2
    xv = 5
    ev = 0.25

    # Act
    ia = IntervalAnalysis(p)
    actual = ia.absolute(context={
        "x": xv
    })


    # Assert
    #assert actual == __bound(expected)