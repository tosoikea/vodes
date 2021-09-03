from abc import ABC, abstractmethod
import logging
import matplotlib.pyplot
#
from vodes.ode.problem import Problem
from numpy.core.function_base import linspace
from sympy import pprint, lambdify

class Solver(ABC):
    # Initialize the solver, this involves the construction of the expression tree
    @abstractmethod
    def init(self):
        pass
    
    # Initialize the solver, this involves the construction of the expression tree
    @abstractmethod
    def parse(self,opt:dict):
        pass

    def __init__(self, problem : Problem, opt = None):
        self.logger = logging.getLogger(__name__)

        # Problem definition including equation and further parameters
        self.problem = problem
        self.solutions = []
        self.precision = None
        # --

        # optional arguments parsed by the implementing solvers
        self.opt = {} if opt is None else opt

        # validate the supplied parameters
        self.parse(opt=opt)
        self.validate()

        # Expression Tree is constructed by the implementing solver
        self.exprs = self.init()

    def validate(self):
        if self.problem is None:
            exit(-1)

        self.problem.validate()

    def calculate(self,precision:int=None):
        from vodes.symbolic.mapper.interop import ExactPymbolicToMathMapper, ExactPymbolicToSympyMapper

        self.symbolic = precision is None
        self.precision = precision

        assert self.symbolic or self.precision > 0

        self.solutions = [
            (
                self.exprs[0][0], 
                ExactPymbolicToSympyMapper()(self.exprs[0][1]) if self.symbolic else ExactPymbolicToMathMapper(precision=precision)(self.exprs[0][1]) 
            )
        ]

        for i in range(1,len(self.exprs)):
            (t,expr) = self.exprs[i]
            (_,y) = self.solutions[-1]
            
            t_approx = ExactPymbolicToSympyMapper()(t) if self.symbolic else ExactPymbolicToMathMapper(precision=precision,context={})(t)
            y_approx = ExactPymbolicToSympyMapper(context={f"y_{i-1}" : y})(expr) if self.symbolic else ExactPymbolicToMathMapper(precision=precision,context={f"y_{i-1}" : y})(expr)

            self.solutions.append(
                (t_approx,y_approx)
            )

    def show(self, exact = None, ticks=100):
        if len(self.solutions) == 0:
            return

        ax = []
        ay = []
        ex = []
        ey = []

        for i in range(0,len(self.solutions)):
            (t,y) = self.solutions[i]
            ax.append(t)
            ay.append(y)

        matplotlib.pyplot.plot(ax, ay)

        if not(exact is None):
            # exact function is no longer dependent of it self
            lexact = lambdify(exact.free_symbols, exact, modules=['numpy'])
            
            # evenly spaced samples of size ticks
            ex = linspace(self.problem.start(),self.problem.end(),num=ticks)
            ey = lexact(ex)
            matplotlib.pyplot.plot(ex, ey)

        matplotlib.pyplot.show()

    # Global error as (absolute) maximum difference between exact and approximated values
    # Obv. not rigorous because of fp
    def error(self, exact, precision:int):
        error = 0

        # exact function is no longer dependent of it self
        lexact = lambdify(exact.free_symbols, exact, modules=['numpy'])
        self.calculate(precision=precision)

        for (t,y) in self.solutions:
            res = lexact(t)
            dif = abs(y - res)
            error = max(dif,error)

        return error

        