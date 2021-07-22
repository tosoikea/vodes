#!/usr/bin/env python

"""Tests for `vodes.error.interval` class."""

from pymbolic.interop.sympy import PymbolicToSympyMapper
import pytest
from random import randrange
from vodes.symbolic.symbols import BoundedExpression, BoundedVariable, DummyVariable
from vodes.symbolic.interval import Interval
from pymbolic import var
from mpmath import mpf,mp
from sympy import lambdify, pprint

def calculate(f,prec):
    t = mp.prec
    try:
        mp.prec = prec
        return f()
    finally:
        mp.prec = t


@pytest.fixture
def analyses():
    from vodes.error.analysis import IntervalAnalysis as IA
    return [
        lambda p : IA(problem=p)
    ]

def __assert_equations(f,bexprs):
    # assume ~ 113 prec is exact (quadruple)
    expected = calculate(f,113)

    # evaluate up to 30 precisions
    ps = [randrange(113) for i in range(randrange(30))]

    for p in ps:
        rounded = calculate(f,p)
        expected_error = abs(expected - rounded)
        print (f'Expected : {expected_error} (prec : {p})')
        
        actual_error = None
        eps = 2 ** (-p)

        for i in range(len(bexprs)):
            bexpr = bexprs[i]

            assert isinstance(bexpr,BoundedExpression)

            if not (bexpr.bound.contains(eps)):
                continue

            f_actual_error = lambdify('eps', PymbolicToSympyMapper()(bexpr.expr))
            actual_error = f_actual_error(eps)

            print (f'Actual : {actual_error} (prec : {p})')
            break

            
        if actual_error is None:
            assert False

        # Actual Error must be positive
        assert actual_error >= 0

        # Actual Error may not be less than the expected error
        assert actual_error >= expected_error


def test_bounds1(analyses):
    # Arrange
    x = var("x")
    c = var("c")
    p = x * c

    cv = randrange(20)
    xv = randrange(20)

    max_prec = 113
    min_prec = 1

    f = lambda: mpf(cv) * mpf(xv)

    for analysis in analyses:
        # Act
        bexprs = analysis(p).absolute(
            context={
                "x" : xv,
                "c" : cv
            },
            max_precision = max_prec,
            min_precision = min_prec
        )

        # Assert
        __assert_equations(
            f = f,
            bexprs=bexprs
        )


def _bounds_test_2(analyses, xv:int, c:int):
    # Arrange
    x = var("x")

    p = (x-1) ** (c) + (1*x/c)

    print(p)
    print(xv)

    max_prec = 113
    min_prec = 1

    f = lambda: (mpf(xv) - 1) ** mpf(c) + (1*mpf(xv) / mpf(c))

    for analysis in analyses:
        # Act
        bexprs = analysis(p).absolute(
            context={
                "x" : xv
            },
            max_precision = max_prec,
            min_precision = min_prec
        )

        # Assert-8.964731982496851e+17 + (1 + e)**3*(779.2 + (1 + e)**4*(3895 + 3896*e)**5))
        __assert_equations(
            f = f,
            bexprs=bexprs
        )


def test_bounds2(analyses):
    c = randrange(start=1,stop=10)
    xv = randrange(200)

    _bounds_test_2(analyses=analyses,xv=xv,c=c)
