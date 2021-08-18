from abc import ABC, abstractmethod
from typing import List

from vodes.symbolic.translations.translation import Translation

# Custom Expression Library
from vodes.symbolic.expressions.bounded import BoundedExpression, Domain
from vodes.symbolic.expressions.interval import Interval

# Expression Library
from pymbolic.primitives import Expression

class NOPTranslation(Translation):
    """Does nothing"""
    def translate(self, expr:BoundedExpression) -> List[BoundedExpression]:
        return [expr]

