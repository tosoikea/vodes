from abc import ABC, abstractmethod
from typing import Type
from mpmath import identification

from sympy import symbols
import functools


def _toFPExpr(f):
    """
    Decorator for converting arguments of a function to FPExpressions before execution.
    """
    @functools.wraps(f)
    def _conversion(*args, **kwargs):
        conv_args = [_convertObjectToFPExpr(a) for a in args]
        return f(*conv_args,**kwargs)
    return _conversion

# Grammar :
# expr = atom | expr o expr
# atom = variable | constant
# 
# with o \in {+,-,/,*}
class FPExpression(ABC):
    def __init__(self):
        if self is FPExpression:
            raise TypeError("FPExpression may not be instantiated")

    @abstractmethod
    def getExpr(self):
        """
        This function is used for obtaining a symbolic representation of the FPExpression, which is in turn stored as tree

        Parameters:
            self (FPExpression): The FPExpression, that is to be converted.
        
        Returns:
            A symbolic expression, that is e.g. a sympy expression. This is mostly for debug purposes.
        """
        pass

    @abstractmethod
    def getIdentification(self):
        """
        This function defines a unique identifier for the given expression. Importantly, two identical expression have to result in the same identifier.

        Parameters:
            self (FPExpression): The FPExpression, that is to be converted to a unique identifier.

        Return:
            An identifier composed of the parameters of the expression.
        """
        pass

    @_toFPExpr
    def __add__(self, other):
        return FPAdd(self,other)

    @_toFPExpr
    def __radd__(self, other):
        return FPAdd(other,self)

    @_toFPExpr
    def __sub__(self, other):
        return FPSub(self,other)

    @_toFPExpr
    def __rsub__(self, other):
        return FPSub(other,self)

    @_toFPExpr
    def __mul__(self, other):
        return FPMul(self,other)

    @_toFPExpr
    def __rmul__(self, other):
        return FPMul(other,self)

    @_toFPExpr
    def __div__(self, other):
        return FPDiv(self,other)

    @_toFPExpr
    def __rdiv__(self, other):
        return FPDiv(other,self)

def _convertObjectToFPExpr(obj):
    """
    Simple conversion from python objects to an equivalent FPExpression, if possible.

    Parameters:
        obj (obj): The object to convert to an FPExpression.
    
    Return:
        A FPExpression, based on the supplied object.
    """
    if isinstance(obj,FPExpression):
        return obj
    elif isinstance(obj, str):
        return FPVar(obj)
    elif isinstance(obj, (int, float, complex)):
        return FPCon(obj)
    else:
        raise TypeError(f'{str(obj)} is not convertable to an FPExpr')


#### Atomic FPExpression

class FPAtom(FPExpression):
    def __init__(self,value):
        if self is FPExpression:
            raise TypeError("FPAtom may not be instantiated")
        self.value = value

    def __str__(self):
        return str(self.value)

    def getIdentification(self):
        return str(self.value).replace(" ", "")

    def __eq__(self, other):
        if isinstance(other, self.__class__):
            return self.value == other.value
        else:
            return False

class FPCon(FPAtom):
    # TODO : limit const
    def __init__(self, const):
        super().__init__(const)

    def getExpr(self):
        return self.value

class FPVar(FPAtom):
    def __init__(self,literal:str):
        super().__init__(literal)

    def getExpr(self):
        return symbols(self.value)
#### ---

#### FPExpression composed of operation

class FPBinaryOperation(FPExpression):
    def __init__(self,left: FPExpression,right: FPExpression):
        if self is FPBinaryOperation:
            raise TypeError("FPBinaryOperation may not be instantiated")
        self._left = left
        self._right = right

    def __str__(self):
        return f'({str(self._left)} {self._op()} {str(self._right)})'

    def getIdentification(self):
        return f'{self._left.getIdentification()}_{self._op()}_{self._right.getIdentification()}'

    @abstractmethod
    def _op(self):
        pass

    def __eq__(self, other):
        if isinstance(other, self.__class__):
            return self._left == other._left and self._right == other._right
        else:
            return False

class FPAdd(FPBinaryOperation):
    def __init__(self,left:FPExpression,right:FPExpression):
        super().__init__(left,right)
        
    def _op(self):
        return "+"

    def getExpr(self):
        return self._left.getExpr() + self._right.getExpr()

class FPMul(FPBinaryOperation):
    def __init__(self,left:FPExpression,right:FPExpression):
        super().__init__(left,right)
        
    def _op(self):
        return "*"

    def getExpr(self):
        return self._left.getExpr() * self._right.getExpr()

class FPSub(FPBinaryOperation):
    def __init__(self,left:FPExpression,right:FPExpression):
        super().__init__(left,right)
        
    def _op(self):
        return "-"

    def getExpr(self):
        return self._left.getExpr() - self._right.getExpr()

class FPDiv(FPBinaryOperation):
    def __init__(self,left:FPExpression,right:FPExpression):
        super().__init__(left,right)
        
    def _op(self):
        return "/"

    def getExpr(self):
        return self._left.getExpr() / self._right.getExpr()