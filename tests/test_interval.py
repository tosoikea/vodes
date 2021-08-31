#!/usr/bin/env python

"""Tests for `vodes.error.interval` class."""
import pytest

from random import randrange

# Custom Expression Library
from vodes.symbolic.expressions.bounded import Domain, BoundedExpression, BoundedVariable, DummyVariable
from vodes.symbolic.expressions.interval import Interval
from vodes.symbolic.expressions.nthroot import NthRoot

# Symbolic Expression Library
from pymbolic import var
from pymbolic.primitives import Quotient, Power

# Utility functions for comparison
from tests.utils import assert_bounded_iv_equations

@pytest.fixture
def evaluators():
    from vodes.symbolic.mapper.comparison_evaluator import ComparisonEvaluator as CE
    
    return [
        lambda c,s : CE(context=c, symbol=s)
    ]

@pytest.fixture
def static_evaluators():
    from vodes.symbolic.mapper.taylor_comparison_evaluator import TaylorComparisonEvaluator as TE
    from vodes.symbolic.mapper.comparison_evaluator import ComparisonEvaluator as CE
    from vodes.symbolic.mapper.scalar_evaluator import ScalarEvaluator as SE

    return [
       SE(context={}, symbol=DummyVariable()),
       CE(context={}, symbol=DummyVariable()),
       TE(context={}, symbol=DummyVariable())
    ]

def assert_static(res,val):
    assert(len(res)==1)
    assert(isinstance(res[0],BoundedExpression))

    assert res[0].expr == val


def test_setup():
    # Arrange
    upper = randrange(5000)
    lower = upper - randrange(5000)

    # Act 
    i = Interval(lower,upper)

    # Assert
    assert i.up == upper
    assert i.low == lower

def test_scalar_add(static_evaluators):
    # Arrange
    upper = randrange(5000)
    lower = upper - randrange(5000)
    a = randrange(5000)

    e = Interval(lower=lower + a, upper=upper + a)

    for eval in static_evaluators:
        # Act
        li = eval(Interval(lower,upper) + a)
        ri = eval(a + Interval(lower,upper))

        # Assert
        assert_static(li, e)
        assert_static(ri, e)

def test_scalar_sub(static_evaluators):
    # Arrange
    upper = randrange(5000)
    lower = upper - randrange(5000)
    a = randrange(5000)

    el = Interval(lower = lower - a, upper = upper - a)
    er = Interval(lower = a - upper, upper = a - lower)

    for eval in static_evaluators:
        # Act
        li = eval(Interval(lower,upper) - a)
        ri = eval(a - Interval(lower,upper))

        # Assert
        assert_static(li, el)
        assert_static(ri, er)

def test_scalar_mul(static_evaluators):
    # Arrange
    upper = randrange(5000)
    lower = upper - randrange(5000)
    a = randrange(5000)

    e = Interval(lower = lower * a, upper = upper * a)

    for eval in static_evaluators:
        # Act
        li = eval(Interval(lower,upper) * a)
        ri = eval(a * Interval(lower,upper))

        # Assert
        assert_static(li, e)
        assert_static(ri, e)

def test_interval_add(static_evaluators):
    # Arrange
    u1 = randrange(5000)
    l1 = u1 - randrange(5000)

    u2 = randrange(5000)
    l2 = u2 - randrange(5000)

    for eval in static_evaluators:
        # Act
        li = eval(Interval(l1,u1) + Interval(l2,u2))
        ri = eval(Interval(l2,u2) + Interval(l1,u1))

        # Assert
        assert_static(li, Interval(lower=l1+l2,upper=u1+u2))
        assert_static(ri, Interval(lower=l1+l2,upper=u1+u2))
    
def test_interval_sub(static_evaluators):
    # Arrange
    u1 = randrange(5000)
    l1 = u1 - randrange(5000)

    u2 = randrange(5000)
    l2 = u2 - randrange(5000)

    for eval in static_evaluators:
        # Act
        li = eval(Interval(l1,u1) - Interval(l2,u2))
        ri = eval(Interval(l2,u2) - Interval(l1,u1))

        # Assert
        assert_static(li, Interval(lower=l1-u2,upper=u1-l2))
        assert_static(ri, Interval(lower=l2-u1,upper=u2-l1))
  
def test_interval_mul(static_evaluators):
    # Arrange
    u1 = randrange(5000)
    l1 = u1 - randrange(5000)

    u2 = randrange(5000)
    l2 = u2 - randrange(5000)

    e = Interval(
        lower=min(u1 * l2, u1 * u2, l1 * l2, l1 * u2),
        upper=max(u1 * l2, u1 * u2, l1 * l2, l1 * u2)
    )

    for eval in static_evaluators:
        # Act
        li = eval(Interval(l1,u1) * Interval(l2,u2))
        ri = eval(Interval(l2,u2) * Interval(l1,u1))

        # Assert
        assert_static(li, e)
        assert_static(ri, e)

def test_interval_pow1(evaluators):
    #Arrange
    u1 = randrange(start=2,stop=22,step=2) #Even exponent
    x = var("x") 
    
    symbol = BoundedVariable(
        x.name,
        boundary=Domain(
            start=-2,
            end=2
        )
    )

    # [x-1,x+1]**u1
    p = Interval(x-1,x+1) ** u1

    # 1. z**u1 has (global) minimum at z == 0
    # 2. f: x - 1 == 0 <=> x = 1 
    #   a) 1 \in [-2,2] (Boundary)
    #   b) 0 <= f([-2,1]), 0 >= f([1,2])
    # 3. f: x + 1 == 0 <=> x = -1
    #   a) -1 \in [-2,2] (Boundary)
    #   b) 0 <= f([-2,-1]), 0 >= f([-1,2])
    # 4. 
    # =>    [-2,-1] min : (x+1)**2 ,   max : (x-1)**2
    #       (-1, 0] min : 0        ,   max : (x-1)**2
    #       (0 , 1] min : 0        ,   max : (x+1)**2
    #       (1 , 2] min : (x-1)**2 ,   max : (x+1)**2

    e = [
        BoundedExpression(
            boundary=Domain(
                start=-2,
                end=-1
            ),
            expression=Interval(
                lower=Power(x+1,u1),
                upper=Power(x-1,u1)
            )
        ),
        BoundedExpression(
            boundary=Domain(
                start=-1,
                end=0,
                left_open=True
            ),
            expression=Interval(
                lower=0,
                upper=Power(x-1,u1)
            )
        ),
        BoundedExpression(
            boundary=Domain(
                start=0,
                end=1,
                left_open=True
            ),
            expression=Interval(
                lower=0,
                upper=Power(x+1,u1)
            )
        ),
        BoundedExpression(
            boundary=Domain(
                start=1,
                left_open=True,
                end=2
            ),
            expression=Interval(
                lower=Power(x-1,u1),
                upper=Power(x+1,u1)
            )
        ),
    ]

    context = {}

    for eval in evaluators:
        #Act
        a = eval(context,symbol)(p)
        print(list(map(str,a)))
        #Assert
        assert_bounded_iv_equations(a,e)

def test_interval_pow2(evaluators):
    #Arrange
    u1 = randrange(1,21,step=2) #Uneven exponent
    x = var("x") 
    
    symbol = BoundedVariable(
        x.name,
        boundary=Domain(
            start=-2,
            end=2,
        )
    )

    p = Interval(x-1,x+1) ** u1

    e = [
        BoundedExpression(
            boundary=Domain(
                start=-2,
                end=2,
            ),
            expression=Interval(
                lower=Power(x-1,u1),
                upper=Power(x+1,u1)
            )
        ),
    ]

    context = {}

    for eval in evaluators:
        #Act
        a = eval(context,symbol)(p)

        #Assert
        assert_bounded_iv_equations(a,e)
        
        
def test_interval_pow3(evaluators):
    #Arrange
    u1 = randrange(1,121)
    u2 = randrange(1,11,step=2) # Uneven exponent
    x = var("x") 
    
    symbol = BoundedVariable(
        x.name,
        boundary=Domain(
            start=0,
            left_open=True,
            end=1,
            right_open=True
        )
    )

    p = Interval(x)**u1 + Interval(-x,x)**u2

    e = [
        BoundedExpression(
            boundary=Domain(
                start=0,
                left_open=True,
                end=1,
                right_open=True
            ),
            expression=Interval(
                lower=x**u1 - x**u2,
                upper=x**u1 + x**u2
            )
        ),
    ]

    context = {}

    for eval in evaluators:
        #Act
        a = eval(context,symbol)(p)

        #Assert
        assert_bounded_iv_equations(a,e)


def test_interval_pow4(evaluators):
    # Arrange
    u1 = 2 #Even exponent
    x = var("x") 
    
    symbol = BoundedVariable(
        x.name,
        boundary=Domain(
            start=-5,
            end=5
        )
    )

    # 1. z**u1 has (global) minimum at z == 0
    # 2. f: x == 0 <=> x = 0 
    #   a) 0 \in [-5,5] (Boundary)
    #   b) 0 <= f([0,5]), 0 >= f([-5,0])
    # 3. f: 0.01x^3-10 == 0 <=> x = 10
    #   a) 10 \not\in [-5,5] (Boundary)
    #   b) 0 >= f([-5,5])
    # 4. 
    # =>    [-5, 0] min : x^2      ,   max : (0.01x^3-10)^2
    #       (0,  5] min : 0        ,   max : (0.01x^3-10)^2

    # [x,0.01*x^3-10]
    lower = Quotient(1,100) * x ** 3 - 10
    i = Interval(lower, x)
    p = i**u1

    e = [
        BoundedExpression(
            boundary=Domain(
                start=-5,
                end=0
            ),
            expression=Interval(
                lower=x**u1,
                upper=lower**u1
            )
        ),
        BoundedExpression(
            boundary=Domain(
                start=0,
                end=5,
                left_open=True
            ),
            expression=Interval(
                lower=0,
                upper=lower**u1
            )
        )
    ]

    context = {}

    for eval in evaluators:
        #Act
        a = eval(context,symbol)(p)

        #Assert
        assert_bounded_iv_equations(a,e)

def test_interval_nthroot1(evaluators):
    # Arrange
    u1 = randrange(2,12,step=2)
    x = var("x") 
    
    symbol = BoundedVariable(
        x.name,
        boundary=Domain(
            start=-100,
            end=100
        )
    )

    # 
    p = NthRoot(Interval(x), u1)**3

    e = [
        BoundedExpression(
            boundary=Domain(
                start=0,
                end=100
            ),
            expression=Interval(
                lower=Power(NthRoot(x,u1),3),
                upper=Power(NthRoot(x,u1),3),
            )
        )
    ]

    context = {}

    for eval in evaluators:
        #Act
        a = eval(context,symbol)(p)

        #Assert
        assert_bounded_iv_equations(a,e)

def test_interval_nthroot2(evaluators):
    # Arrange
    u1 = 2
    x = var("x") 
    
    symbol = BoundedVariable(
        x.name,
        boundary=Domain(
            start=-100,
            end=100
        )
    )

    # 
    p = NthRoot(Interval(x),u1)

    e = [
        BoundedExpression(
            boundary=Domain(
                start=0,
                end=100
            ),
            expression=Interval(
                lower=NthRoot(x,u1),
                upper=NthRoot(x,u1)
            )
        )
    ]

    context = {}

    for eval in evaluators:
        #Act
        a = eval(context,symbol)(p)

        #Assert
        assert_bounded_iv_equations(a,e)

def test_interval_nthroot3(evaluators):
    # Arrange
    u1 = 2
    x = var("x") 
    
    symbol = BoundedVariable(
        x.name,
        boundary=Domain(
            start=0,
            end=100
        )
    )

    # 
    p = NthRoot(Interval(x-1,x+1),u1)

    e = [
        BoundedExpression(
            boundary=Domain(
                start=1,
                end=100
            ),
            expression=Interval(
                lower=NthRoot(x-1,u1),
                upper=NthRoot(x+1,u1)
            )
        )
    ]

    context = {}

    for eval in evaluators:
        #Act
        a = eval(context,symbol)(p)

        #Assert
        assert_bounded_iv_equations(a,e)

def test_interval_div1(evaluators):
    #Arrange
    x = var("x") 
    
    symbol = BoundedVariable(
        x.name,
        boundary=Domain(
            start=-2,
            end=2
        )
    )

    p = 1 / Interval(x-1,x+1)

    e = [
        BoundedExpression(
            boundary=Domain(
                start=-2,
                end=-1,
                right_open=True
            ),
            expression=Interval(
                lower=1/(x+1),
                upper=1/(x-1)
            )
        ),
        BoundedExpression(
            boundary=Domain(
                start=1,
                left_open=True,
                end=2
            ),
            expression=Interval(
                lower=1/(x+1),
                upper=1/(x-1)
            )
        ),
    ]

    context = {}

    for eval in evaluators:
        #Act
        a = eval(context,symbol)(p)
        #Assert
        assert_bounded_iv_equations(a,e)

def test_interval_div2(evaluators):
    #Arrange
    x = var("x") 
    
    symbol = BoundedVariable(
        x.name,
        boundary=Domain(
            start=-1,
            end=1
        )
    )

    p = 1 / Interval(x-1,x+1)

    e = []

    context = {}

    for eval in evaluators:
        #Act
        a = eval(context,symbol)(p)
        #Assert
        assert_bounded_iv_equations(a,e)
        
def test_interval_div3(evaluators):
    #Arrange
    x = var("x") 
    
    symbol = BoundedVariable(
        x.name,
        boundary=Domain(
            start=-5,
            end=5
        )
    )

    p = 1 / Interval(0.5*x**2,x**2)

    e = [
        BoundedExpression(
            boundary=Domain(
                start=-5,
                end=0,
                right_open=True
            ),
            expression=Interval(
                lower=1/(x**2),
                upper=1/(0.5*x**2)
            )
        ),
        BoundedExpression(
            boundary=Domain(
                start=0,
                left_open=True,
                end=5
            ),
            expression=Interval(
                lower=1/(x**2),
                upper=1/(0.5*x**2)
            )
        )
    ]

    context = {}

    for eval in evaluators:
        #Act
        a = eval(context,symbol)(p)
        #Assert
        assert_bounded_iv_equations(a,e)

def test_interval_div4(evaluators):
    #Arrange
    x = var("x") 
    
    symbol = BoundedVariable(
        x.name,
        boundary=Domain(
            start=-5,
            end=5
        )
    )

    p = 1 / Interval(0.5*x**2,x+20)

    e = [
        BoundedExpression(
            boundary=Domain(
                start=-5,
                end=0,
                right_open=True
            ),
            expression=Interval(
                lower=1/(x+20),
                upper=1/(0.5*x**2)
            )
        ),
        BoundedExpression(
            boundary=Domain(
                start=0,
                left_open=True,
                end=5
            ),
            expression=Interval(
                lower=1/(x+20),
                upper=1/(0.5*x**2)
            )
        )
    ]

    context = {}

    for eval in evaluators:
        #Act
        a = eval(context,symbol)(p)
        #Assert
        assert_bounded_iv_equations(a,e)

def test_interval_div5(evaluators):
    import sympy as sp
    #Arrange
    x = var("x") 
    
    symbol = BoundedVariable(
        x.name,
        boundary=Domain(
            start=-1,
            end=1,
            right_open=True
        )
    )

    expr = x**4 - x**2 + Power(4,-1)
    p = 1 / Interval(expr)

    e = [
        BoundedExpression(
            boundary=Domain(
                start=-1,
                end=-sp.sqrt(2)/2,
                right_open=True
            ),
            expression=Interval(
                lower=1/(expr),
                upper=1/(expr)
            )
        ),
        BoundedExpression(
            boundary=Domain(
                start=-sp.sqrt(2)/2,
                left_open=True,
                end=sp.sqrt(2)/2,
                right_open=True
            ),
            expression=Interval(
                lower=1/(expr),
                upper=1/(expr)
            )
        ),
        BoundedExpression(
            boundary=Domain(
                start=sp.sqrt(2)/2,
                left_open=True,
                end=1,
                right_open=True
            ),
            expression=Interval(
                lower=1/(expr),
                upper=1/(expr)
            )
        )
    ]

    context = {}

    for eval in evaluators:
        #Act
        a = eval(context,symbol)(p)
        #Assert
        assert_bounded_iv_equations(a,e)

def test_interval_div6(evaluators):
    import sympy as sp
    #Arrange
    x = var("x") 
    
    symbol = BoundedVariable(
        x.name,
        boundary=Domain(
            start=-100,
            end=100,
            right_open=True
        )
    )

    p1 = 1 / Interval(x)
    p2 = Interval(x)**(-1)

    context = {}

    for eval in evaluators:
        #Act
        a1 = eval(context,symbol)(p1)
        a2 = eval(context,symbol)(p2)

        #Assert
        assert_bounded_iv_equations(a1,a2)



