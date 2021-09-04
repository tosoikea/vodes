import logging
from abc import ABC, abstractmethod


# Custom Expressions
from vodes.symbolic.expressions.absolute import Abs
from vodes.symbolic.expressions.bounded import  BoundedExpression, MachineError, SmallestMachineNumber
from vodes.symbolic.expressions.interval import Interval

# Custom Mapper
from vodes.error.mapper import IntervalMapper
from vodes.symbolic.mapper.binary_mapper import BinaryMapper as BM
from vodes.symbolic.mapper.interop import ExactPymbolicToSympyMapper, ExactSympyToPymbolicMapper
from vodes.symbolic.mapper.taylor_comparison_evaluator import TaylorComparisonEvaluator as IE

# Symbolic Expression
from pymbolic.primitives import Expression, Power

# Mapper

class Analysis(ABC):
    """Superclass for the concise implementations of roundoff error analysis for a symbolic problem.

    Args:
        problem (Expression): The problem, that is to be analyzed for its roundoff errors.

    Attributes:
        _problem (Expression): The supplied problem stored in binary format. Every Node therefore has a maximum of two children.
    """
    def __init__(self, problem:Expression,opt:dict):
        self._opt = {} if opt is None else opt
        self._logger = logging.getLogger(__name__)
        self._problem = BM()(problem)
        self._absolute = None
        self._precision = None

        self.parse(self._opt)

    @abstractmethod
    def parse(self, opt:dict):
        pass

    @abstractmethod
    def absolute(self, context:dict, min_precision:int, max_precision:int, min_exponent:int, max_exponent:int) -> list:
        """Calculate the absolute error for the earlier supplied problem.

        Args:
            context (dict): Defines the substitution of symbolic variables within the problem for the evaluation. Importantly, all free variables contained within the initial problem HAVE to be substituted. Otherwise, the execution terminates.
            min_precision (int): Defines the minimum precision that should be used for the analysis. 
            max_precision (int): Defines the maximum precision that should be used for the analysis.
            min_exponent (int): Minimum bits used for the exponent.
            max_exponent (int): Maximum bits used for the exponent.
            
        Raises:
            UnknownVariableError

        Returns:
            absolute: Either a list of symbolic expressions, if the machine error ("e") was not substituted or the absolute error as a constant. This is a rigerous approximation (greater, but not less).
        """
        pass

    def show(self):
        """Visualize the calculated absolute error."""
        from vodes.error.utils import show, AnalysisSolution

        if self._absolute is None:
            raise ValueError("No absolute error calculated")

        show(
            solutions=[
                AnalysisSolution(bexprs=self._absolute,name=str(self.__class__))
            ],
            min_prec=self._precision[0],
            max_prec=self._precision[1]
        )

class IntervalAnalysis(Analysis):
    """Class implementing roundoff error analysis on the basis of interval analysis.

    Args:
        problem (Expression): The problem, that is to be analyzed for its roundoff errors.

    Attributes:
        _problem (Expression): The supplied problem stored in binary format. Every Node therefore has a maximum of two children.
        __expr (Expression): The problem with appropriate error terms injected.
    """
    def __init__(self, problem: Expression, opt:dict=None):
        super().__init__(problem, opt)

        # Mapper exposes all created noise variables
        self._mapper = IntervalMapper

    def parse(self, opt:dict):
        # No optional parameters
        return

    def absolute(self, context:dict, min_precision:int, max_precision:int, min_exponent:int, max_exponent:int) -> list:
        # (1) Create Floating Point Error Model
        self.__mapper_instance = IntervalMapper(
            context=context,
            min_precision=min_precision,
            max_precision=max_precision,
            min_exponent=min_exponent,
            max_exponent=max_exponent
        )
        _expr = self.__mapper_instance(self._problem)

        err = Abs(self._problem - _expr)

        self._precision = (min_precision,max_precision)
        self._absolute = IE(
            context=context,
            symbol=MachineError(
                min_precision=min_precision,
                max_precision=max_precision
                )
            )(err)

        self._logger.info(f'Evaluated Error : {list(map(str,self._absolute))}')

        return self._absolute

class TaylorAnalysis(Analysis):
    """Class implementing roundoff error analysis on the basis of taylor expansions
    
    Args:
        problem (Expression): The problem, that is to be analyzed for its roundoff errors.
        symbolic (bool): Determines if the remainder of the taylor expressions should be calculated using epsilon as a free variable.

    Attributes:
        _problem (Expression): The supplied problem stored in binary format. Every Node therefore has a maximum of two children.
        __expr (Expression): The problem with appropriate error terms injected.
    """
    def __init__(self, problem: Expression, opt: dict=None):
        super().__init__(problem, opt)

    def parse(self, opt:dict):
        self.n = opt.setdefault("n", 1)
        self._logger.info(f"Requested taylor order is {self.n}")

    def absolute(self, context:dict, min_precision:int, max_precision:int, min_exponent:int, max_exponent:int) -> list:
        from vodes.symbolic.mapper.extended_substitution_mapper import substitute
        from vodes.symbolic.mapper.comparison_evaluator import evaluate
        from vodes.error.mapper import TaylorMapper, expand_taylor_terms
        import copy

        self.__mapper_instance = TaylorMapper(
            n=self.n,
            context=copy.deepcopy(context),
            min_precision=min_precision,
            max_precision=max_precision,
            min_exponent=min_exponent,
            max_exponent=max_exponent
        )
        
        (_, error, _) = expand_taylor_terms(
                self.__mapper_instance(self._problem),
                min_prec=min_precision,
                max_prec=max_precision,
                abs=True
            )

        self._absolute = evaluate(Abs(substitute(
            error,
            context
        )),min_precision,max_precision)

        self._logger.info(f"Determined absolute error {list(map(str,self._absolute))}")
        return self._absolute