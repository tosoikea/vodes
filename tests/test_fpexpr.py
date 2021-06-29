#!/usr/bin/env python

"""Tests for `vodes.error.interval` class."""

from random import randrange
import pytest
from vodes.error.fpexpr import *

def test_eq1():
    # Arrange
    a = randrange(5000)

    # Act
    c1 = FPCon(a)
    c2 = FPCon(a)

    # Assert
    assert c1 == c2

def test_eq2():
    # Arrange
    a = randrange(5000)

    # Act
    c1 = FPVar(str(a))
    c2 = FPVar(str(a))

    # Assert
    assert c1 == c2
    
def test_eq2():
    # Arrange
    a = randrange(5000)
    c1 = FPVar(a)
    c2 = FPCon(a)


    # Act
    aa1 = FPAdd(c1,c2)
    aa2 = FPAdd(c1,c2)
    ab1 = FPAdd(c2,c1)
    ab2 = FPAdd(c2,c1)

    sa1 = FPSub(c1,c2)
    sa2 = FPSub(c1,c2)
    sb1 = FPSub(c2,c1)
    sb2 = FPSub(c2,c1)

    ma1 = FPMul(c1,c2)
    ma2 = FPMul(c1,c2)
    mb1 = FPMul(c2,c1)
    mb2 = FPMul(c2,c1)

    da1 = FPDiv(c1,c2)
    da2 = FPDiv(c1,c2)
    db1 = FPDiv(c2,c1)
    db2 = FPDiv(c2,c1)

    # Assert
    assert aa1 == aa2
    assert ab1 == ab2
    assert sa1 == sa2
    assert sb1 == sb2
    assert ma1 == ma2
    assert mb1 == mb2
    assert da1 == da2
    assert db1 == db2

def test_conversion_a1():
    # Arrange
    a = FPVar("a")
    b = randrange(5000)

    # Act
    r = a + b

    # Assert
    assert r == FPAdd(a,FPCon(b))

def test_conversion_a2():
    # Arrange
    a = FPVar("a")
    b = randrange(5000)

    # Act
    r = b + a

    # Assert
    assert r == FPAdd(FPCon(b),a)

def test_conversion_s1():
    # Arrange
    a = FPVar("a")
    b = randrange(5000)

    # Act
    r = a - b

    # Assert
    assert r == FPSub(a,FPCon(b))

def test_conversion_s2():
    # Arrange
    a = FPVar("a")
    b = randrange(5000)

    # Act
    r = b - a

    # Assert
    assert r == FPSub(FPCon(b),a)

def test_conversion_m1():
    # Arrange
    a = FPVar("a")
    b = randrange(5000)

    # Act
    r = a * b

    # Assert
    assert r == FPMul(a,FPCon(b))


def test_conversion_m2():
    # Arrange
    a = FPVar("a")
    b = randrange(5000)

    # Act
    r = b * a

    # Assert
    assert r == FPMul(FPCon(b),a)

def test_conversion_c1():
    # Arrange
    y = FPVar("y")
    h = FPVar("h")
    c = FPCon(3)

    # Act
    f = y + h * c * y

    # Assert
    assert f == FPAdd(y,FPMul(FPMul(h,c),y))