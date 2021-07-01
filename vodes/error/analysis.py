from abc import ABC, abstractmethod

from pymbolic.mapper.coefficient import CoefficientCollector
from vodes.error.mapper import IntervalMapper
from vodes.symbolic.binary_mapper import BinaryMapper as BM
from pymbolic.primitives import Expression
from pymbolic.functions import fabs
from vodes.error.ia_evaluation import MachineEvaluator, BoundEvaluator


from pymbolic.mapper.constant_folder import ConstantFoldingMapper

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
        err = self._problem - self.__expr

        noised = MachineEvaluator(context=context)(err)
        print("N : ")
        print(noised)

        #bounded = BoundEvaluator(context=context)(noised)
        #print("B : ")
        #print(bounded)
        print("F")
        print(CoefficientCollector()(noised))
        print("C")
        print(ConstantFoldingMapper()(noised))

        return MachineEvaluator(context=context)(err)