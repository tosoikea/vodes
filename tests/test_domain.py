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
    
    