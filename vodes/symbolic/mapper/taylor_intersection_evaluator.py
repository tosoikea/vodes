from typing import List, Tuple
from vodes.symbolic.translations.to_taylor import ToTaylor
from vodes.symbolic.properties.is_polynomial import IsPolynomial
from vodes.symbolic.analysis import Analysis

# Assumption library
from vodes.symbolic.translations.to_scalar import ToScalar
from vodes.symbolic.properties.is_scalar import IsScalar
from vodes.symbolic.assumption import Assumption

# Custom Expression library
from vodes.symbolic.expressions.interval import Interval
from vodes.symbolic.expressions.bounded import BoundedVariable, BoundedExpression, Domain
from vodes.symbolic.utils import le,ge,gt,lt,minimum,maximum

# Custom Mappers
from vodes.symbolic.mapper.comparison_evaluator import ComparisonEvaluator

# Expression Library
from pymbolic.primitives import Expression, Quotient, Variable, Power

class TaylorIntersectionEvaluator(ComparisonEvaluator):
    """Class for determining the exact boundaries of intervals on the basis of function analysis.
    TODO : Describe in depth, firstly in thesis."""

    def __init__(self, context: dict, symbol: BoundedVariable):
        super().__init__(context=context, symbol=symbol)
        
        self._assumptions["_minimum"] = [
                Assumption(
                    property=IsPolynomial(),
                    translation=ToTaylor(n=2)
                )
            ]

        self._assumptions["_maximum"] = [
                Assumption(
                    property=IsPolynomial(),
                    translation=ToTaylor(n=2)
                )
            ]
    