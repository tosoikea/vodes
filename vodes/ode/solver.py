from abc import ABC, abstractmethod
import logging
import matplotlib.pyplot
#
from vodes.ode.problem import Problem
from numpy.core.function_base import linspace
from sympy import pprint, lambdify

class Solver(ABC):
    # The symbolic encoding of the solver. MUST be initiated in the child classes
    # All free variables, that are substituted later on, must be included in the returned tuple
    # e.g. (y + h * 3 * y,[ y, h ])
    @abstractmethod
    def _sym_step(self):
        pass

    # Return parameters for lambdafied function
    @abstractmethod
    def _get_params(self, y ,t):
        pass

    def __init__(self, problem : Problem, opt = None):
        self.logger = logging.getLogger(__name__)

        # Precision of floating point calculation
        # TODO : Proper usage within calculation
        self.precision = 24

        # Problem definition including equation and further parameters
        self.problem = problem
        self.solution = []

        # steps continously appended or initially constructed by the implementing solver
        self.steps = []
        # --

        # optional arguments parsed by the implementing solvers
        self.opt = {} if opt is None else opt

        # validate the supplied parameters
        self.validate()

        # Initialize solver
        step, symbols = self._sym_step()

        self.__lambda_step = lambdify(symbols, step)
        # --

    def validate(self):
        if self.problem is None:
            exit(-1)

        self.problem.validate()

    ## way more efficient than symbolic evaluation
    def _eval(self, y, t):
        return self.__lambda_step(*self._get_params(y ,t))

    def store(self, t, y):
        self.logger.debug('Calculated %f for %d',y,t)
        self.solution.append((t,y))

    def calculate(self):
        if len(self.solution) != 0:
            self.solution = []
            self.steps = []

        # initial values
        self.solution.append((self.problem.start(), self.problem.initial()))
        self.logger.debug('Starting calculate with %f at %d',self.solution[0][1], self.solution[0][0])

        while len(self.steps) > 0:
            t = self.steps.pop(0)
            y = self.solution[-1][1]

            res = self._eval(self.solution[-1][0],y)
            self.store(t, res)

    def show(self, exact = None, ticks=100):
        if len(self.solution) == 0:
            return

        ax = []
        ays = [[]]
        ex = []
        ey = []

        for i in range(0,len(self.solution)):
            ax.append(self.solution[i][0])

            if self.problem.matrix:
                for j in range(0, len(self.solution[i][1])):
                    if len(ays) <= j:
                        ays.append([])
                    ays[j].append(self.solution[i][1][j])
            else:
                ays[0].append(self.solution[i][1])


        if not(exact is None):
            # exact function is no longer dependent of it self
            lexact = lambdify(exact.free_symbols, exact, modules=['numpy'])
            
            # evenly spaced samples of size ticks
            ex = linspace(self.problem.start(),self.problem.end(),num=ticks)
            ey = lexact(ex)
            matplotlib.pyplot.plot(ex, ey)

        for ay in ays:
            matplotlib.pyplot.plot(ax, ay)
        matplotlib.pyplot.show()

    # Global error as (absolute) maximum difference between exact and approximated values
    def error(self, exact):
        error = 0

        # exact function is no longer dependent of it self
        lexact = lambdify(exact.free_symbols, exact, modules=['numpy'])

        for i in range(0, len(self.solution)):
            res = lexact(self.solution[i][0])
            dif = abs(self.solution[i][1] - res)
            error = max(dif,error)

        return error
        
        