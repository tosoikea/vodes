from abc import ABC, abstractmethod

from numpy.core.function_base import linspace
from vodes.symbolic.absolute import Abs
from vodes.symbolic.maximum import Max

from pymbolic.mapper.evaluator import EvaluationMapper
from pymbolic.interop.sympy import PymbolicToSympyMapper
from vodes.symbolic.symbols import BoundedExpression, MachineError
from vodes.symbolic.interval import ExactIntersectionEvaluator

from vodes.error.mapper import IntervalMapper
from vodes.symbolic.binary_mapper import BinaryMapper as BM
from pymbolic.primitives import Expression
from pymbolic.functions import fabs
from sympy import lambdify
from matplotlib import pyplot
from math import ceil
import logging

class Analysis(ABC):
    """Superclass for the concise implementations of roundoff error analysis for a symbolic problem.

    Args:
        problem (Expression): The problem, that is to be analyzed for its roundoff errors.

    Attributes:
        _problem (Expression): The supplied problem stored in binary format. Every Node therefore has a maximum of two children.
    """
    def __init__(self, problem:Expression):
        self._logger = logging.getLogger(__name__)
        self._problem = BM()(problem)
        self._absolute = None

    @abstractmethod
    def absolute(self, context:dict, min_precision:int, max_precision:int) -> list:
        """Calculate the absolute error for the earlier supplied problem.

        Args:
            context (dict): Defines the substitution of symbolic variables within the problem for the evaluation. Importantly, all free variables contained within the initial problem HAVE to be substituted. Otherwise, the execution terminates.
            min_precision (int): Defines the minimum precision that should be used for the analysis. 
            max_precision (int): Defines the maximum precision that should be used for the analysis.
            
        Raises:
            UnknownVariableError

        Returns:
            absolute: Either a list of symbolic expressions, if the machine error ("e") was not substituted or the absolute error as a constant. This is a rigerous approximation (greater, but not less).
        """
        pass

    def show(self, ticks=100, end=1):
        if self._absolute is None:
            raise ValueError("No absolute error calculated")

        min_start = self._absolute[0].bound.lower.value
        max_end = self._absolute[-1].bound.upper.value

        for bexpr in self._absolute:
            # TODO : handle open intervals
            start = EvaluationMapper()(bexpr.bound.lower.value)
            end = min(EvaluationMapper()(bexpr.bound.upper.value),end)

            if end <= start:
                continue
            
            bticks = ceil(ticks * (end - start) / (max_end - min_start))

            expr_err = PymbolicToSympyMapper()(bexpr.expr)
            f_err = lambdify(expr_err.free_symbols, expr_err)

            # TODO : absolute error may have no free symbols, if machine symbol was passed
            xs = linspace(start, end, num=bticks)
            ys = f_err(xs)

            pyplot.plot(xs,ys)

        pyplot.yscale('log',base=2)
        pyplot.xscale('log',base=2)
        pyplot.show()


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

        self._logger.info(f'Machine Expression : {self.__expr}')

    def absolute(self, context: dict,min_precision:int,max_precision:int):
        if min_precision <= 0 or max_precision < min_precision:
            raise ValueError(f"The supplied precision values {min_precision} and {max_precision} are invalid")

        err = Max(Abs(self._problem - self.__expr))
        self._logger.info(f'Error : {err}')

        self._absolute = ExactIntersectionEvaluator(
            context=context,
            symbol=MachineError(
                min_precision=min_precision,
                max_precision=max_precision
                )
            )(err)

        self._logger.info(f'Evaluated Error : {list(map(str,self._absolute))}')

        return self._absolute