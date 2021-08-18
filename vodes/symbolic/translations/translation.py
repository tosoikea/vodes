from abc import ABC, abstractmethod
from typing import List

# Custom Expression Library
from vodes.symbolic.expressions.bounded import BoundedExpression, Domain

class Translation(ABC):
    """Superclass for translating bounded expressions with the goal of obtaining certain properties."""

    @abstractmethod
    def translate(self, expr:BoundedExpression) -> List[BoundedExpression]:
        """Translates the given expression bounded by the domain."""
        pass
        
