from abc import ABC, abstractmethod
import logging

class Solver(ABC):
    def __init__(self, problem, y0, steps):
        self.problem = problem
        self.y0 = y0
        self.steps = steps
        self.solution = []

        self.logger = logging.getLogger(__name__)

        # validate the supplied parameters
        self.validate()

    @abstractmethod
    def validate(self):
        pass

    @abstractmethod
    def step(self, y, t):
        pass

    def calculate(self):
        while len(self.steps) > 0:
            t = self.steps.pop(0)
            y = self.y0 if len(self.solution) == 0 else self.solution[-1][1]
            res = self.step(y,t)

            self.logger.info('Calculated %f for %d',res,t)

            self.solution.append((t,res))