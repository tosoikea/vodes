from abc import ABC, abstractmethod

from pymbolic.interop.sympy import PymbolicToSympyMapper
from vodes.symbolic.symbols import MachineError
from vodes.symbolic.interval import ExactIntersectionEvaluator

from vodes.error.mapper import IntervalMapper
from vodes.symbolic.binary_mapper import BinaryMapper as BM
from pymbolic.primitives import Expression
from pymbolic.functions import fabs

class Analysis(ABC):
    def __init__(self, problem:Expression):
        self._problem = BM()(problem)

    @abstractmethod
    def absolute(self, context:dict):
        pass


class IntervalAnalysis(Analysis):
    def __init__(self, problem: Expression):
        super().__init__(problem)
        self.__expr = IntervalMapper()(self._problem)

    def absolute(self, context: dict):
        print(self.__expr)
        #err = fabs(self._problem - self.__expr)
        err = self._problem - self.__expr

        noised = ExactIntersectionEvaluator(context=context,symbol=MachineError())(err)
        print("N : ")
        print(noised)
        return noised