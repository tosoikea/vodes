from abc import ABC, abstractmethod
from typing import List

# Custom Expression Library
from vodes.symbolic.expressions.bounded import BoundedExpression

class Property(ABC):
    """Superclass for the properties of symbolic expressions, that can be validated."""

    @abstractmethod
    def verify(self, expr:BoundedExpression) -> bool:
        """Verifies the property for the given expression bounded by the domain."""
        pass
        
