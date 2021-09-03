from vodes.ode.solver import Solver
from pymbolic.primitives import Variable
from abc import ABC, abstractmethod

class RK(Solver,ABC):
    def __init__(self, problem, opt = None):
        super().__init__(problem, opt = opt)

    def parse(self,opt:dict):
        self.dt = opt.setdefault("dt", 1)
        self.logger.debug("Fixed step size is %f", self.dt)

    def validate(self):
        from vodes.symbolic.utils import gt
        super().validate()

        assert gt(self.dt,0)

    def init(self):
        from vodes.symbolic.utils import le
        steps = []

        while True:
            step = self.problem.start
            if len(steps) > 0:
                step = steps[-1]

            tn = step + self.dt
            if le(tn, self.problem.end):
                steps.append(tn)
            else:
                break
        
        iterations = len(steps)
        self.logger.info(f"Determined {len(steps)} steps from {self.problem.start} to {self.problem.end} with step size {self.dt}")
        
        t = self.problem.start
        exprs = [
            (t,self.problem.initial)
        ]

        for i in range(iterations):
            y = self._step(Variable(f'y_{i}'),t)
            t = t + self.dt 
            exprs.append(
                (t,y)
            )

        self.logger.info(f'Finished setting up expression.')
        return exprs

    @abstractmethod
    def _step(self, y, t):
        pass

class RK1(RK):
    def __init__(self, problem, opt = None):
        super().__init__(problem, opt = opt)

    def _step(self, y, t):
        # y_n + h * f(y_n,t_n)
        return y + self.dt * self.problem.get_expression(y, t)

class RK2(RK):
    def _step(self, y, t):
        from vodes.symbolic.mapper.extended_substitution_mapper import substitute

        k1 = Variable('k1')
        k2 = Variable('k2')
        
        # y_n + h / 2 * [ f(y_n,t_n) + f(t_n + h, y_n + h * f(y_n + t_n)) ]
        step = y + self.dt / 2 * (k1 + k2)

        step = substitute(step, {k2 : self.problem.get_expression(y + self.dt * k1, t + self.dt)})
        step = substitute(step, {k1 : self.problem.get_expression(y, t)})

        return step

class RK4(RK):
    def _kutta_step(self, problem):
        k1, k2, k3, k4 = self._get_sym(["k1", "k2", "k3", "k4"])

        step = self.y + self.h / 6 * (k1 + 2 * k2 + 2 * k3 + k4)

        step = step.subs(k4, problem.get_expression(self.y + self.h * k3, self.t + self.h))
        step = step.subs(k3, problem.get_expression(self.y + self.h / 2 * k2, self.t + self.h / 2))
        step = step.subs(k2, problem.get_expression(self.y + self.h / 2 * k1, self.t + self.h / 2))
        step = step.subs(k1, problem.get_expression(self.y, self.t))

        return step
