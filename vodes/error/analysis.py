import logging
from abc import ABC, abstractmethod


# Custom Expressions
from vodes.symbolic.expressions.absolute import Abs
from vodes.symbolic.expressions.bounded import  BoundedExpression, MachineError, SmallestMachineNumber
from vodes.symbolic.expressions.interval import Interval

# Custom Mapper
from vodes.error.mapper import AffineMapper, IntervalMapper
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
    def __init__(self, problem:Expression):
        self._logger = logging.getLogger(__name__)
        self._problem = BM()(problem)
        self._absolute = None
        self._precision = None

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
    def __init__(self, problem: Expression,symbolic:bool=True):
        from vodes.symbolic.mapper.scalar_evaluator import ScalarEvaluator
        super().__init__(problem)

        # Mapper exposes all created noise variables
        self._mapper = IntervalMapper
        self._evaluator = IE if symbolic else ScalarEvaluator

    def absolute(self, context:dict, min_precision:int, max_precision:int, min_exponent:int, max_exponent:int) -> list:
        # (1) Create Floating Point Error Model
        self.__mapper_instance = self._mapper(
            context=context,
            min_precision=min_precision,
            max_precision=max_precision,
            min_exponent=min_exponent,
            max_exponent=max_exponent
        )
        _expr = self.__mapper_instance(self._problem)

        err = Abs(self._problem - _expr)
        self._logger.info(f'Error : {err}')

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
    def __init__(self, problem: Expression, symbolic:bool = False):
        super().__init__(problem)

        # Mapper exposes all created noise variables
        self._mapper = AffineMapper
        self._symbolic = symbolic

    def _create_error_terms(self,context,min_precision:int,max_precision:int,min_exponent:int,max_exponent:int):
        from sympy import diff
        from functools import reduce
        from vodes.symbolic.analysis import Analysis as FuncAnalysis
        from vodes.symbolic.mapper.interop import ExactSympyToPymbolicMapper as STP

        # (1) Create Floating Point Error Model
        self.__mapper_instance = self._mapper(
            context=context,
            min_precision=min_precision,
            max_precision=max_precision,
            min_exponent=min_exponent,
            max_exponent=max_exponent
        )
        _expr = self.__mapper_instance(self._problem)
        _analysis = FuncAnalysis(_expr)

        self._logger.info(f'Floating Point Model : {_expr}')

        # (2) Create Taylor expansions
        # We do not multiply with (x-a), as it is approximated later on

        ##
        # Absolute Error := | o(f(x)) - f(x) |
        # Affine form with noise variables to model rounded operations
        #                <= | f(x,e) - f(x) |
        # Taylor relative to e at e = (0,...,0)
        # T_1 : first order term
        # R_1 : first order error term
        #                <= | f(x,0) + T_1(f,x,0) + R_1(f,x,p) - f(x) |
        # f(x,0) = f(x)
        #                 = | T_1(f,x,0) + R_1(f,x,p) |
        ##
        _t1 = _analysis.taylor(n=1,no_tail=True)

        ##
        # Taylor’s Inequality :
        # M >= | f^(n+1)(x) | for all x \in (- \eps, + \eps) => | R_n(x) | <= M / (n+1)! * | x - a |^(n+1)
        #
        # R_1(x) = (x - a)^T / 2 * Hf(a)(x - a)
        #
        # d^2 f / de^2  = \sum_{i=0,j=0} (\partial^2 f) / (\partial e_i \partial e_j) (x,e)
        # https://mathinsight.org/taylors_theorem_multivariable_introduction
        #
        ##
        _r1 = _analysis.remainder(n=1)

        return (STP()(_t1),STP()(_r1))

    def absolute(self, context:dict, min_precision:int, max_precision:int, min_exponent:int, max_exponent:int) -> list:
        from vodes.symbolic.mapper.extended_substitution_mapper import substitute
        from vodes.symbolic.utils import merge

        # (1) Taylor Expansion Terms
        (_t1,_r1) = self._create_error_terms(context=context,min_precision=min_precision,max_precision=max_precision,min_exponent=min_exponent,max_exponent=max_exponent)

        print(f'T1 : {_t1}')
        print(f'R1 : {_r1}')

        # (2) Substitutions for Remainder evaluation

        # (A) Maximum Relative Error
        eps = MachineError(min_precision=min_precision,max_precision=max_precision)

        # (B) Maximum Absolute Error for subnormal calculations
        sigma = SmallestMachineNumber(min_exponent=min_exponent,max_exponent=max_exponent)

        subs = {
            "sigma": Interval(sigma.bound.start,sigma.bound.end)
        }

        for key in self.__mapper_instance.exact:
            subs[key] = Interval(-eps, eps) if self._symbolic else Interval(-eps.bound.end, eps.bound.end)

        # (1) Remainder, may or may not be symbolic
        # TODO : Also expose exponent, we use it as a constant (after substitution)
        eps_r1 = IE(context=context,symbol=eps)(
            Abs(substitute(_r1,subs))
        )

        # (2) Linear Expression, does not contain symbolic expressions
        eps_t1 = IE(context=context,symbol=eps)(
            Abs(_t1)
        )  

        # (3) Final Expression
        self._precision = (min_precision,max_precision)

        self._absolute  = [
            BoundedExpression(
                # TODO : Sigma as free variable
                expression=eps * et1 + eps**2 * er1,
                boundary=boundary
            ) for ((et1,er1),boundary) in merge(eps_t1,eps_r1)
        ]
        self.__mapper_instance = None

        self._logger.info(f'Evaluated Error : {list(map(str,self._absolute))}')
        return self._absolute
