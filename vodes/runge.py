from vodes.solver import Solver
from abc import ABC, abstractmethod

class RK(Solver,ABC):
    def validate(self):
        super().validate()

        self.dt = self.opt.setdefault("dt", 1)
        self.logger.debug("Fixed step size is %f", self.dt)

        self._next_step(self.interval[0])

    # Prepares next step
    def _next_step(self, t):
        tn = self.dt + t
        if tn <= self.interval[1]:
            self.steps.append(tn)

    def store(self, t, y):
        super().store(t, y)
        self._next_step(t)

class RK1(RK):
    def step(self, t, y):
        k1 = self._eval(t,y)
        return y + self.dt * k1 

class RK2(RK):
    def step(self, t, y):
        k1 = self._eval(t,y)
        k2 = self._eval(t + self.dt, y + self.dt * k1)
        return y + self.dt / 2 * (k1 + k2)
