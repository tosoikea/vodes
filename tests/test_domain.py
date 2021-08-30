#!/usr/bin/env python

"""Tests for `vodes.error.interval` class."""
import pytest

from random import randrange

# Custom Expression Library
from vodes.symbolic.expressions.bounded import Domain

def test_difference1():
    # (1,2] - [1,2) = [2,2]

    # Arrange
    l = Domain(1,2,left_open=True,right_open=False)
    r = Domain(1,2,left_open=False,right_open=True)

    # Act
    res = l.difference([r])

    # Assert
    expected = [
        Domain(
            2,2,False,False
        )
    ]

    assert expected == res
    
def test_difference2():
    # [1,2] - [-1,1] = (1,2]

    # Arrange
    l = Domain(1,2,left_open=False,right_open=False)
    r = Domain(-1,1,left_open=False,right_open=False)

    # Act
    res = l.difference([r])

    # Assert
    expected = [
        Domain(
            1,2,True,False
        )
    ]

    assert expected == res

def test_difference3():
    # [2,4] - (2,4) = [2,2],[4,4]

    # Arrange
    l = Domain(2,4,left_open=False,right_open=False)
    r = Domain(2,4,left_open=True,right_open=True)

    # Act
    res = l.difference([r])

    # Assert
    expected = [
        Domain(
            2,2,False,False
        ),
        Domain(
            4,4,False,False
        )
    ]

    assert expected == res

def test_difference4():
    # [0,4] - [1,4) = [0,1),[4,4]

    # Arrange
    l = Domain(0,4,left_open=False,right_open=False)
    r = Domain(1,4,left_open=False,right_open=True)

    # Act
    res = l.difference([r])

    # Assert
    expected = [
        Domain(
            0,1,False,True
        ),
        Domain(
            4,4,False,False
        )
    ]

    assert expected == res

def test_difference5():
    # [-1,2] - [-100,100) = []

    # Arrange
    l = Domain(-1,2,left_open=False,right_open=False)
    r = Domain(-100,100,left_open=False,right_open=True)

    # Act
    res = l.difference([r])

    # Assert
    expected = []

    assert expected == res

def test_union1():
    # (-1,1] u [1,1] = [(-1,1]]

    # Arrange
    l = Domain(-1,1,left_open=True,right_open=False)
    r = Domain(1,1,left_open=False,right_open=False)

    # Act
    r1 = l.union(r)
    r2 = r.union(l)

    # Assert
    expected = [l]

    assert expected == r1
    assert expected == r2

def test_union2():
    # (-1,1] u (1,3] = [(-1,3]]

    # Arrange
    l = Domain(-1,1,left_open=True,right_open=False)
    r = Domain(1,3,left_open=True,right_open=False)

    # Act
    r1 = l.union(r)
    r2 = r.union(l)

    # Assert
    e1 = [Domain(-1,3,True,False)]

    assert e1 == r1
    assert e1 == r2

def test_union3():
    # [-1,1) u (1,3] = [[-1,1),(1,3]]

    # Arrange
    l = Domain(-1,1,left_open=False,right_open=True)
    r = Domain(1,3,left_open=True,right_open=False)

    # Act
    r1 = l.union(r)
    r2 = r.union(l)

    # Assert
    e1 = [l,r]
    e2 = [r,l]

    assert e1 == r1
    assert e2 == r2

def test_union4():
    # [-1,0] u [1,3] = [[-1,0],[1,3]]

    # Arrange
    l = Domain(-1,0,left_open=False,right_open=False)
    r = Domain(1,3,left_open=False,right_open=False)

    # Act
    r1 = l.union(r)
    r2 = r.union(l)

    # Assert
    e1 = [l,r]
    e2 = [r,l]

    assert e1 == r1
    assert e2 == r2

def test_union5():
    # [-1,0) u [0,3) = [[-1,3)]

    # Arrange
    l = Domain(-1,0,left_open=False,right_open=True)
    r = Domain(0,3,left_open=False,right_open=True)

    # Act
    r1 = l.union(r)
    r2 = r.union(l)

    # Assert
    e1 = [Domain(-1,3,False,True)]

    assert e1 == r1
    assert e1 == r2

def test_union6():
    # [-1,20) u [10,15] = [[-1,20)]

    # Arrange
    l = Domain(-1,20,left_open=False,right_open=True)
    r = Domain(10,15,left_open=False,right_open=False)

    # Act
    r1 = l.union(r)
    r2 = r.union(l)

    # Assert
    e1 = [l]

    assert e1 == r1
    assert e1 == r2

def test_union7():
    # [-1,20) u [10,25] = [[-1,25]]

    # Arrange
    l = Domain(-1,20,left_open=False,right_open=True)
    r = Domain(10,25,left_open=False,right_open=False)

    # Act
    r1 = l.union(r)
    r2 = r.union(l)

    # Assert
    e1 = [Domain(-1,25,False,False)]

    assert e1 == r1
    assert e1 == r2
    
    
    