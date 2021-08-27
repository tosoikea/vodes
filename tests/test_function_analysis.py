#!/usr/bin/env python
import pytest
"""Tests for `vodes.symbolic.analysis`."""

from random import randrange

from vodes.symbolic.analysis import Analysis, AnalysisConfig
from vodes.symbolic.expressions.interval import Interval
from vodes.symbolic.mapper.interop import ExactPymbolicToSympyMapper

# Custom Expression Library
from vodes.symbolic.expressions.bounded import Domain, BoundedExpression

# Symbolic Expression Library
from pymbolic import var


def __assert_inclusion(expr,bexprs):
    from vodes.symbolic.expressions.infinity import NegativeInfinity, Infinity
    from sympy import solveset
    from sympy.core.numbers import NegativeInfinity as NegInf, Infinity as Inf
    from sympy.sets.sets import EmptySet, Interval as SymInterval

    for bexpr in bexprs:
        assert isinstance(bexpr, BoundedExpression)
        assert isinstance(bexpr.expr, Interval)

        boundary = SymInterval(
            NegInf() if isinstance(bexpr.bound.start,NegativeInfinity) else bexpr.bound.start,
            Inf() if isinstance(bexpr.bound.end,Infinity) else bexpr.bound.end,
            bexpr.boundary.left_open,
            bexpr.boundary.right_open
        )

        print(f'Evaluating inclusion of {expr} in {bexpr.expr} within {boundary}')

        sym_expr = ExactPymbolicToSympyMapper()(expr)
        sym_low = ExactPymbolicToSympyMapper()(bexpr.expr.low)
        sym_up = ExactPymbolicToSympyMapper()(bexpr.expr.up)

        # (1) expr >= iv.low?
        res = solveset(sym_expr < sym_low,list(sym_expr.free_symbols)[0],boundary)
        assert isinstance(res,EmptySet)

        # (2) iv.up >= expr?
        res = solveset(sym_expr > sym_up,list(sym_expr.free_symbols)[0],boundary)
        assert isinstance(res,EmptySet)

def test_inclusion_const1():
    """Simple inclusion of constant expression"""
    # Arrange
    x = 0
    an = Analysis(x)

    # Act
    bexprs = an.taylor(n=randrange(1,9))

    # Assert
    assert len(bexprs) == 1
    assert isinstance(bexprs[0],BoundedExpression)
    assert x == bexprs[0].expr.low
    assert x == bexprs[0].expr.up

def test_inclusion_const2():
    """Simple inclusion of constant expression that is randomly generated"""
    # Arrange
    x = randrange(-1000,1000)
    an = Analysis(x)

    # Act
    bexprs = an.taylor(n=randrange(1,9))

    # Assert
    assert len(bexprs) == 1
    assert isinstance(bexprs[0],BoundedExpression)
    assert x == bexprs[0].expr.low
    assert x == bexprs[0].expr.up

def test_inclusion_linear1():
    """Simple inclusion of linear expression"""
    from tests.utils import __assert_equivalent

    # Arrange
    x = var('x')
    a = 2
    b = 200
    p = a * x + b

    an = Analysis(p)

    # Act
    bexprs = an.taylor(n=randrange(1,9))

    # Assert
    assert len(bexprs) == 1
    assert isinstance(bexprs[0],BoundedExpression)
    __assert_equivalent(p,bexprs[0].expr.low)
    __assert_equivalent(p,bexprs[0].expr.up)
    
def test_inclusion_linear2():
    """Simple inclusion of linear expression that is randomly generated"""
    from tests.utils import __assert_equivalent

    # Arrange
    x = var('x')
    a = randrange(-3203,3203)
    b = randrange(-3203,3203)
    p = a * x + b

    an = Analysis(p)

    # Act
    bexprs = an.taylor(n=randrange(1,9))

    # Assert
    assert len(bexprs) == 1
    assert isinstance(bexprs[0],BoundedExpression)
    __assert_equivalent(p,bexprs[0].expr.low)
    __assert_equivalent(p,bexprs[0].expr.up)

def test_inclusion_quadratic1():
    """Inclusion of quadratic expression"""
    # Arrange
    x = var('x')
    a1 = 2
    a2 = 4
    b = 200
    p = b + a1 * x + a2 * x**2
    d = None

    an = Analysis(p,config=AnalysisConfig(d=d))

    # Act
    res = [
        an.taylor(n=1),
        an.taylor(n=2),
        an.taylor(n=3)
    ]

    # Assert
    for bexprs in res:
        __assert_inclusion(p,bexprs)

def test_inclusion_quadratic2():
    """Inclusion of quadratic expression that is randomly generated"""
    # Arrange
    x = var('x')
    a1 = randrange(-2000,3000)
    a2 = randrange(-2000,3000)
    b = randrange(-2000,3000)
    p = b + a1 * x + a2 * x**2
    d = None

    an = Analysis(p,config=AnalysisConfig(d=d))

    # Act
    res = [
        an.taylor(n=1),
        an.taylor(n=2),
        an.taylor(n=3)
    ]

    # Assert
    for bexprs in res:
        __assert_inclusion(p,bexprs)

def test_inclusion_quadratic3():
    """Inclusion of quadratic expression that lead to failures in earlier runs"""
    # Arrange
    x = var('x')
    a0 = -1655
    a1 = -1188
    a2 = -1288
    p = a0 + a1 * x + a2 * x**2
    d = None

    an = Analysis(p,config=AnalysisConfig(d=d))

    # Act
    res = [
        an.taylor(n=1),
        an.taylor(n=2),
        an.taylor(n=3)
    ]

    # Assert
    for bexprs in res:
        __assert_inclusion(p,bexprs)

def test_inclusion_cubic1():
    """Inclusion of cubic expression"""
    # Arrange
    x = var('x')
    a0 = 33
    a1 = 1
    a2 = 3
    a3 = 9
    p = a0 + a1 * x + a2 * x**2 + a3 * x**3
    d = None

    an = Analysis(p,config=AnalysisConfig(d=d))

    # Act & Assert
    try:
        an.taylor(n=1)
        assert False
    except:
        assert True

def test_inclusion_cubic2():
    """Inclusion of cubic expression"""
    # Arrange
    x = var('x')
    a0 = 33
    a1 = 1
    a2 = 3
    a3 = 9
    p = a0 + a1 * x + a2 * x**2 + a3 * x**3
    d = Domain(0,1)

    an = Analysis(p,config=AnalysisConfig(d=d))

    # Act
    res = [
        an.taylor(n=1),
        an.taylor(n=2),
        an.taylor(n=3)
    ]

    # Assert
    for bexprs in res:
        __assert_inclusion(p,bexprs)

def test_inclusion_cubic3():
    """Inclusion of cubic expression that was randomly generated"""
    # Arrange
    x = var('x')
    a0 = randrange(-2000,3000)
    a1 = randrange(-2000,3000)
    a2 = randrange(-2000,3000)
    a3 = randrange(-2000,3000)
    p = a0 + a1 * x + a2 * x**2 + a3 * x**3
    d = Domain(randrange(-2,2),randrange(2,6))

    an = Analysis(p,config=AnalysisConfig(d=d))

    # Act
    res = [
        an.taylor(n=1),
        an.taylor(n=2),
        an.taylor(n=3),
        an.taylor(n=4)
    ]

    # Assert
    for bexprs in res:
        __assert_inclusion(p,bexprs)

def test_inclusion_cubic4():
    """Inclusion of cubic expression that lead to a failure in previous runs"""
    # Arrange
    x = var('x')
    a0 = 1043
    a1 = -1110
    a2 = 1654
    a3 = -1536
    p = a0 + a1 * x + a2 * x**2 + a3 * x**3
    d = Domain(randrange(-2,2),randrange(2,6))

    an = Analysis(p,config=AnalysisConfig(d=d))

    # Act
    res = [
        an.taylor(n=1),
        an.taylor(n=2),
        an.taylor(n=3),
        an.taylor(n=4)
    ]

    # Assert
    for bexprs in res:
        __assert_inclusion(p,bexprs)

def test_inclusion_cubic5():
    """Inclusion of cubic expression that lead to a failure in previous runs"""
    # Arrange
    x = var('x')
    a0 = 2039
    a1 = -860
    a2 = -1754
    a3 = 1609
    p = a0 + a1 * x + a2 * x**2 + a3 * x**3
    d = Domain(randrange(-2,2),randrange(2,6))

    an = Analysis(p,config=AnalysisConfig(d=d))

    # Act
    res = [
        an.taylor(n=1),
        an.taylor(n=2),
        an.taylor(n=3),
        an.taylor(n=4)
    ]

    # Assert
    for bexprs in res:
        __assert_inclusion(p,bexprs)

def test_inclusion_quartic1():
    """Inclusion of quartic expression"""
    # Arrange
    x = var('x')
    a0 = 33
    a1 = 1
    a2 = 3
    a3 = 9
    a4 = 27
    p = a0 + a1 * x + a2 * x**2 + a3 * x**3 + a4 * x**4
    d = None

    an = Analysis(p,config=AnalysisConfig(d=d))

    # Act & Assert
    try:
        an.taylor(n=1)
        assert False
    except:
        assert True

    try:
        an.taylor(n=2)
        assert False
    except:
        assert True

def test_inclusion_quartic2():
    """Inclusion of quartic expression"""
    # Arrange
    x = var('x')
    a0 = 33
    a1 = 1
    a2 = 3
    a3 = 9
    a4 = 27
    p = a0 + a1 * x + a2 * x**2 + a3 * x**3 + a4 * x**4
    d = Domain(0,1)

    an = Analysis(p,config=AnalysisConfig(d=d))

    # Act
    res = [
        an.taylor(n=1),
        an.taylor(n=2),
        an.taylor(n=3),
        an.taylor(n=4),
        an.taylor(n=5)
    ]

    # Assert
    for bexprs in res:
        __assert_inclusion(p,bexprs)

def test_inclusion_quartic3():
    """Inclusion of quartic expression that was randomly generated"""
    # Arrange
    x = var('x')
    a0 = randrange(-2000,3000)
    a1 = randrange(-2000,3000)
    a2 = randrange(-2000,3000)
    a3 = randrange(-2000,3000)
    a4 = randrange(-2000,3000)
    p = a0 + a1 * x + a2 * x**2 + a3 * x**3 + a4 * x**4
    d = Domain(randrange(-2,2),randrange(2,6))

    an = Analysis(p,config=AnalysisConfig(d=d))

    # Act
    res = [
        an.taylor(n=1),
        an.taylor(n=2),
        an.taylor(n=3),
        an.taylor(n=4),
        an.taylor(n=5)
    ]

    # Assert
    for bexprs in res:
        __assert_inclusion(p,bexprs)

def test_inclusion_quartic4():
    """Inclusion of quartic expression that lead to a failure in previous runs"""
    # Arrange
    x = var('x')
    a0 = -643
    a1 = 343
    a2 = 609
    a3 = 2449
    a4 = 891
    p = a0 + a1 * x + a2 * x**2 + a3 * x**3 + a4 * x**4
    d = Domain(randrange(-2,2),randrange(2,6))

    an = Analysis(p,config=AnalysisConfig(d=d))

    # Act
    res = [
        an.taylor(n=1),
        an.taylor(n=2),
        an.taylor(n=3),
        an.taylor(n=4),
        an.taylor(n=5)
    ]

    # Assert
    for bexprs in res:
        __assert_inclusion(p,bexprs)