from abc import ABC, abstractmethod
import logging
import matplotlib.pyplot

class Solver(ABC):
    def __init__(self, problem, y0, interval, opt = None):
        self.problem = problem
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

    @abstractmethod
    def step(self, t, y):
        pass

    def store(self, t, y):
        self.logger.debug('Calculated %f for %d',y,t)
        self.solution.append((t,y))

    def calculate(self):
        if len(self.solution) != 0:
            return

        # initial value
        self.solution.append((self.interval[0], self.y0))
        self.logger.debug('Starting calculate with %f at %d',self.solution[0][1], self.solution[0][0])

        while len(self.steps) > 0:
            t = self.steps.pop(0)
            y = self.solution[-1][1]

            res = self.step(self.solution[-1][0],y)
            self.store(t, res)

    def show(self, exact = None):
        if len(self.solution) == 0:
            return

        x = []
        ay = []
        ey = []

        for i in range(0,len(self.solution)):
            x.append(self.solution[i][0])
            ay.append(self.solution[i][1])

            if exact is None:
                continue

            ey.append(exact(x[-1]))

        if len(ey) == 0:
            matplotlib.pyplot.plot(x, ay)
        else:
            matplotlib.pyplot.plot(x, ay, x, ey)
        matplotlib.pyplot.show()


        
        