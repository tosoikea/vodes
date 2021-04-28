from abc import ABC, abstractmethod
import logging
import matplotlib.pyplot

#
from numpy.core.function_base import linspace
from sympy.core.symbol import symbols
from sympy.utilities.lambdify import lambdify

class Solver(ABC):
    def __init__(self, problem, y0, interval, opt = None, symbols = None):
        self.problem = problem

        # problems are dependent on time and itself
        # TODO : Move Problem into separate class
        self.__internal_symbols = [] if symbols is None else symbols
        self.__internal_problem = lambdify(self.__internal_symbols, self.problem)

        self.y0 = y0
        self.interval = interval
        self.solution = []

        # optional arguments parsed by the implementing solvers
        self.opt = {} if opt is None else opt

        # steps continously appended or initially constructed by the implementing solver
        self.steps = []

        self.logger = logging.getLogger(__name__)


        # validate the supplied parameters
        self.validate()

    def validate(self):
        if len(self.interval) != 2:
            exit(-1)

        if len(self.__internal_symbols) > 2:
            exit(-1)

    def _eval(self, t, y):
        if len(self.__internal_symbols) == 0:
            return self.__internal_problem()
        elif len(self.__internal_symbols) == 1:
            return self.__internal_problem(y)
        else:
            return self.__internal_problem(y, t)

    @abstractmethod
    def step(self, t, y):
        pass

    def store(self, t, y):
        self.logger.debug('Calculated %f for %d',y,t)
        self.solution.append((t,y))

    def calculate(self):
        if len(self.solution) != 0:
            self.solution = []
            self.steps = []
            self.validate()

        # initial value
        self.solution.append((self.interval[0], self.y0))
        self.logger.debug('Starting calculate with %f at %d',self.solution[0][1], self.solution[0][0])

        while len(self.steps) > 0:
            t = self.steps.pop(0)
            y = self.solution[-1][1]

            res = self.step(self.solution[-1][0],y)
            self.store(t, res)

    def show(self, exact = None, ticks=100):
        if len(self.solution) == 0:
            return

        ax = []
        ay = []
        ex = []
        ey = []

        for i in range(0,len(self.solution)):
            ax.append(self.solution[i][0])
            ay.append(self.solution[i][1])

        if not(exact is None):
            # exact function is no longer dependent of it self
            lexact = lambdify(exact.free_symbols, exact, modules=['numpy'])
            
            # evenly spaced samples of size ticks
            ex = linspace(self.interval[0],self.interval[1],num=ticks)
            ey = lexact(ex)
            matplotlib.pyplot.plot(ex, ey)

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
        
        