#!/usr/bin/env python

"""Tests for `vodes.error.interval` class."""

from vodes.symbolic.mapper.interop import ExactPymbolicToSympyMapper
import pytest
from random import randrange
from vodes.symbolic.expressions.bounded import BoundedExpression
from pymbolic import var
from mpmath import mpf,mp
from sympy import Pow


LOWEST_PRECISION = 11
HIGHEST_PRECISION = 53
MIN_EXP = 8
MAX_EXP = 8

# assume ~ 113 prec is exact (quadruple)
EXACT_PRECISION = 113


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
    from vodes.error.analysis import TaylorAnalysis as TA
    return [
        lambda p : IA(problem=p),
        #lambda p : TA(problem=p)
    ]

def __assert_equations(f,bexprs):
    expected = calculate(f,EXACT_PRECISION)

    # evaluate up to 30 precisions for half to double precision
    ps = [randrange(LOWEST_PRECISION,HIGHEST_PRECISION) for i in range(randrange(30))]

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
            
            actual_error = ExactPymbolicToSympyMapper()(bexpr.expr).subs("eps", Pow(2,-p))
            print (f'Actual : {actual_error} (prec : {p})')
            break

        if actual_error is None:
            assert False

        # Actual Error must be positive
        assert actual_error >= 0

        # Actual Error may not be less than the expected error
        assert actual_error >= expected_error


def test_inequality1(analyses):
    # Arrange
    x = var("x")
    c = var("c")
    p = x * c

    cv = randrange(20)
    xv = randrange(20)

    f = lambda: mpf(cv) * mpf(xv)

    for analysis in analyses:
        # Act
        bexprs = analysis(p).absolute(
            context={
                "x" : xv,
                "c" : cv
            },
            max_precision = HIGHEST_PRECISION,
            min_precision = LOWEST_PRECISION,
            min_exponent = MIN_EXP,
            max_exponent = MAX_EXP
        )

        # Assert
        __assert_equations(
            f = f,
            bexprs=bexprs
        )


def _inequality2(analyses, xv:int, c:int):
    # Arrange
    x = var("x")

    p = (x-1) ** (c) + (1*x/c)

    f = lambda: (mpf(xv) - 1) ** mpf(c) + (1*mpf(xv) / mpf(c))

    for analysis in analyses:
        # Act
        bexprs = analysis(p).absolute(
            context={
                "x" : xv
            },
            max_precision = HIGHEST_PRECISION,
            min_precision = LOWEST_PRECISION,
            min_exponent = MIN_EXP,
            max_exponent = MAX_EXP
        )

        # Assert
        __assert_equations(
            f = f,
            bexprs=bexprs
        )

def test_inequality2a(analyses):
    c = randrange(start=1,stop=10)
    xv = randrange(200)

    _inequality2(analyses=analyses,xv=xv,c=c)
    
def test_inequality2b(analyses):
    """This test case was problematic in earlier iterations"""
    c = 5
    xv = 36

    _inequality2(analyses=analyses,xv=xv,c=c)

def test_inequality2c(analyses):
    """This test case was problematic in earlier iterations (NoConvergence)"""
    c = 2
    xv = 60

    _inequality2(analyses=analyses,xv=xv,c=c)

def test_inequality3(analyses):
    from vodes.symbolic.expressions.nthroot import NthRoot
    from vodes.symbolic.mapper.extended_evaluation_mapper import evaluate
    from pymbolic.primitives import Quotient
    from mpmath import sqrt

    # Arrange
    x = var('x')
    sym_xv = Quotient(1,10*6)
    xv = evaluate(sym_xv)

    p = NthRoot(1+x,2) - 1


    f = lambda: sqrt(mpf(xv) + mpf(1)) - mpf(1)

    for analysis in analyses:
        # Act
        bexprs = analysis(p).absolute(
            context={
                "x" : sym_xv
            },
            max_precision = HIGHEST_PRECISION,
            min_precision = LOWEST_PRECISION,
            min_exponent = MIN_EXP,
            max_exponent = MAX_EXP
        )

        # Assert
        __assert_equations(
            f = f,
            bexprs=bexprs
        )

# Does not work because of https://github.com/sympy/sympy/issues/20097
# Wait for version 1.9 and replay test
def test_inequality5(analyses):
    from vodes.symbolic.expressions.nthroot import NthRoot
    from mpmath import sqrt

    return

    # Arrange
    x = var('x')
    xv = 420

    # sqrt(x+1) - sqrt(x)
    p = NthRoot(x+1,n=2) - NthRoot(x,n=2)
    f = lambda: sqrt(mpf(xv) + mpf(1)) - sqrt(mpf(xv))

    for analysis in analyses:
        # Act
        bexprs = analysis(p).absolute(
            context={
                "x" : xv
            },
            max_precision = HIGHEST_PRECISION,
            min_precision = LOWEST_PRECISION
        )

        # Assert
        __assert_equations(
            f = f,
            bexprs=bexprs
        )

# Does not work because of https://github.com/sympy/sympy/issues/20097
# Wait for version 1.9 and replay test
def test_inequality6(analyses):
    from vodes.symbolic.expressions.nthroot import NthRoot
    from mpmath import sqrt
    
    return

    # Arrange
    x = var('x')
    xv = 420

    # 1 / (sqrt(x+1) + sqrt(x))
    p = 1 / (NthRoot(x+1,n=2) + NthRoot(x,n=2))
    f = lambda: mpf(1) / (sqrt(mpf(xv) + mpf(1)) + sqrt(mpf(xv)))

    for analysis in analyses:
        # Act
        bexprs = analysis(p).absolute(
            context={
                "x" : xv
            },
            max_precision = HIGHEST_PRECISION,
            min_precision = LOWEST_PRECISION
        )

        # Assert
        __assert_equations(
            f = f,
            bexprs=bexprs
        )

