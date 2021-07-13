from abc import ABC, abstractmethod
from vodes.symbolic.absolute import Abs

from pymbolic.interop.sympy import PymbolicToSympyMapper
from vodes.symbolic.symbols import MachineError
from vodes.symbolic.interval import ExactIntersectionEvaluator

from vodes.error.mapper import IntervalMapper
from vodes.symbolic.binary_mapper import BinaryMapper as BM
from pymbolic.primitives import Expression
from pymbolic.functions import fabs

class Analysis(ABC):
    """Superclass for the concise implementations of roundoff error analysis for a symbolic problem.

    Args:
        problem (Expression): The problem, that is to be analyzed for its roundoff errors.

    Attributes:
        _problem (Expression): The supplied problem stored in binary format. Every Node therefore has a maximum of two children.
    """
    def __init__(self, problem:Expression):
        self._problem = BM()(problem)

    @abstractmethod
    def absolute(self, context:dict) -> list:
        """Calculate the absolute error for the earlier supplied problem.

        Args:
            context (dict): Defines the substitution of symbolic variables within the problem for the evaluation. Importantly, all free variables contained within the initial problem HAVE to be substituted. Otherwise, the execution terminates.

        Raises:
            UnknownVariableError

        Returns:
            absolute: Either a list of symbolic expressions, if the machine error ("e") was not substituted or the absolute error as a constant. This is a rigerous approximation (greater, but not less).
        """
        pass


class IntervalAnalysis(Analysis):
    """Class implementing roundoff error analysis on the basis of interval analysis.

    Args:
        problem (Expression): The problem, that is to be analyzed for its roundoff errors.

    Attributes:
        _problem (Expression): The supplied problem stored in binary format. Every Node therefore has a maximum of two children.
        __expr (Expression): The problem with appropriate error terms injected.
    """
    def __init__(self, problem: Expression):
        super().__init__(problem)
        self.__expr = IntervalMapper()(self._problem)

        print("Machine Expression :")
        print(self.__expr)

    def absolute(self, context: dict):
        err = Abs(self._problem - self.__expr)
        noised = ExactIntersectionEvaluator(context=context,symbol=MachineError())(err)
        print("N : ")
        print(PymbolicToSympyMapper()(noised[0].expr.up))
        return noised