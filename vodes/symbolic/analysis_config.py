from typing import Type

# Custom Expression
from vodes.symbolic.expressions.bounded import Domain, DummyVariable
from vodes.symbolic.expressions.infinity import Infinity, NegativeInfinity

# Custom Mapper
from vodes.symbolic.mapper.symbolic_interval_evaluator import SymbolicIntervalEvaluator

class AnalysisConfig:
    @property
    def dom(self):
        """Get or set the lower boundary of the interval"""
        return self.__domain

    @property
    def denominator_limit(self):
        """Get or set the lower boundary of the interval"""
        return self.__limit

    def is_multivariant(self):
        """Get the support for multivariant expression"""
        return self.__evaluator.is_multivariant()

    def evaluator(self,context,symbols:list):
        """Get the symbolic interval evaluator"""
        if len(symbols) > 1:
            if not self.is_multivariant():
                raise ValueError(f"Not supporting multivariant expressions with the evaluator {self.__evaluator}")
            else:
                return self.__evaluator(context=context,symbols=symbols)
        elif len(symbols) == 1:
            return self.__evaluator(context=context,symbol=symbols[0])
        else:
            return self.evaluator(context=context,symbols=[DummyVariable()])

    def __init__(self, limit:int=None,d:Domain=None,evaluator:Type[SymbolicIntervalEvaluator]=None):
        from vodes.symbolic.mapper.comparison_evaluator import ComparisonEvaluator

        # TODO : Complex support
        self.__domain = Domain(NegativeInfinity(),Infinity(),left_open=True,right_open=True) if d is None else d
        self.__limit = 10**3 if limit is None else limit

        # We use an interval evaluator, to determine proper symbolic interval constructions.
        # However, they are limited to only support single variant expressions!
        if evaluator is None:
            evaluator = ComparisonEvaluator

        self.__evaluator = evaluator

        assert self.__limit > 0