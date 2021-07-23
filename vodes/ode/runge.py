from vodes.ode.solver import Solver
from sympy import symbols, MatrixSymbol, shape
from abc import ABC, abstractmethod

class RK(Solver,ABC):
    def __init__(self, problem, opt = None):
        self.y = None
        self.t = None
        self.h = None
        super().__init__(problem, opt = opt)

    def validate(self):
        super().validate()

        self.dt = self.opt.setdefault("dt", 1)
        self.logger.debug("Fixed step size is %f", self.dt)

        self._next_step(self.problem.start())

    # Prepares next step
    def _next_step(self, t):
        tn = self.dt + t
        if tn <= self.problem.end():
            self.steps.append(tn)

    def store(self, t, y):
        super().store(t, y)
        self._next_step(t)

    @abstractmethod
    def _kutta_step(self, problem):
        pass

    def _get_params(self, t, y):
        return (y, t, self.dt)

    def _sym_step(self):
        self.y = self._get_sym(["y"])
        self.t, self.h = symbols('t h')
        return (self._kutta_step(self.problem).doit(), [self.y, self.t, self.h])

    def _get_sym(self, names : list):
        if self.problem.matrix:
            if self.y is None:
                (r,_) = shape(self.problem.get_expression(*symbols('y t')))
                c = 1
            else:
                r = self.y.rows
                c = self.y.cols
            res = []
            for name in names:
                res.append(MatrixSymbol(name, r, c))
        else:
            res = symbols(names)

        if len(res) == 1:
            return res[0]
        else:
            return tuple(res)

class RK1(RK):
    def _kutta_step(self, problem):
        # y_n + h * f(y_n,t_n)
        return self.y + self.h * problem.get_expression(self.y, self.t)

class RK2(RK):
    def _kutta_step(self, problem):
        k1, k2 = self._get_sym(["k1", "k2"])
        
        # y_n + h / 2 * [ f(y_n,t_n) + f(t_n + h, y_n + h * f(y_n + t_n)) ]
        step = self.y + self.h / 2 * (k1 + k2)

        step = step.subs(k2, problem.get_expression(self.y + self.h * k1, self.t + self.h))
        step = step.subs(k1, problem.get_expression(self.y, self.t))

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
