from vodes.solver import Solver

def _check_equidistant(steps):
    if len(steps) <= 1:
        exit(-1)
    else:
        dt = abs(steps[0] - steps[1])
        for i in range(2,len(steps)):
            if abs(steps[i - 1] - steps[i]) != dt:
                exit(-1)
        return dt

class RK1(Solver):
    def validate(self):
        super().validate()
        self.dt = _check_equidistant(self.steps)

    def step(self, y, t):
        k1 = self.problem(t,y)
        return y + self.dt * k1

class RK2(Solver):
    def validate(self):
        super().validate()
        self.dt = _check_equidistant(self.steps)

    def step(self, y, t):
        k1 = self.problem(t,y)
        k2 = self.problem(t + self.dt, y + self.dt * k1)
        return y + self.dt / 2 * (k1 + k2)
