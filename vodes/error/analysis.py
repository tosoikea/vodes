import logging
from abc import ABC, abstractmethod


# Custom Expressions
from vodes.symbolic.expressions.absolute import Abs
from vodes.symbolic.expressions.bounded import  MachineError
from vodes.symbolic.expressions.interval import Interval

# Custom Mapper
from vodes.error.mapper import AffineMapper, IntervalMapper
from vodes.symbolic.mapper.binary_mapper import BinaryMapper as BM
from vodes.symbolic.mapper.interop import ExactPymbolicToSympyMapper
from vodes.symbolic.mapper.taylor_intersection_evaluator import TaylorIntersectionEvaluator as IE

# Symbolic Expression
from pymbolic.primitives import Expression, Power

# Mapper
from pymbolic.mapper.substitutor import SubstitutionMapper, make_subst_func
from pymbolic.mapper.evaluator import EvaluationMapper
from pymbolic.interop.sympy import SympyToPymbolicMapper

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

    def show(self, ticks=100):
        from matplotlib import pyplot
        from sympy import Float
        from numpy import linspace
        from math import ceil

        #from numpy.core.function_base import linspace

        if self._absolute is None:
            raise ValueError("No absolute error calculated")

        max_v = None

        for bexpr in self._absolute:
            # TODO : handle open intervals
            from vodes.symbolic.mapper.extended_evaluation_mapper import evaluate
            start = evaluate(bexpr.bound.start)
            end = evaluate(bexpr.bound.end)
            
            # TODO : Implement scaling based on sub-domain length
            bticks = ticks

            expr_err = ExactPymbolicToSympyMapper()(bexpr.expr)

            xs = linspace(start, end, num=bticks)

            if len(expr_err.free_symbols) == 0:
                y = expr_err.evalf(512)
                ys = [y for x in xs]
            elif len(expr_err.free_symbols) == 1:
                sym = expr_err.free_symbols.pop()
                ys = [expr_err.subs(sym,Float(x,512)).evalf(512) for x in xs]
            else:
                raise ValueError("Ecountered too many free variables.")

            pyplot.plot(xs,ys)

            max_v = max(ys) if max_v is None else max(max(ys),max_v)

        pyplot.xscale('log',base=2)
        pyplot.grid(True)
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

        err = Abs(self._problem - self.__expr)
        self._logger.info(f'Error : {err}')

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

    Attributes:
        _problem (Expression): The supplied problem stored in binary format. Every Node therefore has a maximum of two children.
        __expr (Expression): The problem with appropriate error terms injected.
    """
    def __init__(self, problem: Expression):
        super().__init__(problem)
        
        self.__create_error_term()


    def __create_error_term(self):
        from sympy import diff
        from functools import reduce

        # Mapper exposes all created noise variables
        self.__mapper = AffineMapper()

        self.__expr = self.__mapper(self._problem)
        sym_expr = ExactPymbolicToSympyMapper()(self.__expr)

        self._logger.info(f'Affine Expression : {sym_expr}')
        self._logger.debug(f'Mappings : {self.__mapper.ids}')

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

        # We do not multiply with noise (* noise), as we approximate it later, by MachineError
        self.__t1 = SympyToPymbolicMapper()(
            reduce(lambda x,y: x + y, [
                    diff(sym_expr, noise.name).subs(self.__mapper.exact) for noise in self.__mapper.noises
                ]
            )
        )
        self._logger.info(f'T1 : {self.__t1}')

        ##
        # Taylor’s Inequality :
        # M >= | f^(n+1)(x) | for all x \in (- \eps, + \eps) => | R_n(x) | <= M / (n+1)! * | x - a |^(n+1)
        #
        # R_1(x) = (x - a)^T / 2 * Hf(a)(x - a)
        #
        # d^2 f / de^2  = \sum_{i=0,j=0} (\partial^2 f) / (\partial e_i \partial e_j) (x,e)
        # https://mathinsight.org/taylors_theorem_multivariable_introduction
        #
        # Idea : Exact Remainder, Compute Intervall for noises, upper bound is M
        ##

        # We do not multiply with noises (* e_i * e_j), as we approximate it later, by MachineError
        self.__r1 = SympyToPymbolicMapper()(
            reduce(lambda x,y: x + y, [
                    diff(sym_expr,e_i.name,e_j.name) for e_i in self.__mapper.noises for e_j in self.__mapper.noises
                ]
            )
        )
        self._logger.info(f'R1 : {self.__r1}')



    def absolute(self, context: dict,min_precision:int,max_precision:int):
        eps = MachineError(min_precision=min_precision,max_precision=max_precision)

        subs = {}
        for key in self.__mapper.exact:
            subs[key] = Interval(-eps, eps)

        subs_mapper = SubstitutionMapper(subst_func=make_subst_func(subs))
        eps_r1 = subs_mapper(self.__r1)

        
        ##
        # The approximative error term consists no machine error
        ##
        t_t1 = EvaluationMapper(context=context)(self.__t1)
        err = eps * Abs(t_t1) + Power(eps,2) * Abs(eps_r1)
        self._logger.info(f'Approximative Error : {err}')

        self._absolute = IE(
            context=context,
            symbol=eps
        )(err)

        self._logger.info(f'Evaluated Error : {list(map(str,self._absolute))}')

        return self._absolute
