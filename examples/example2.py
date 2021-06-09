#
import sys
sys.path.append('../vodes')
#

from vodes.ode.runge import RK1, RK2, RK4
from vodes.ode.problem import Problem
import matplotlib.pyplot

from sympy import symbols
from sympy.functions import exp

iv = 1
c = 3.0

# define function
t = symbols('t')
y = symbols('y')
f = c * y

problem = Problem(f, iv, (0,1))

# exact solution
e = exp(c * t) * iv

# visualize
steps = []
error_r1 = []
error_r2 = []
error_r3 = []

for i in range(0,6):
    step = 10 ** -i
    print(step)
    lr1 = RK1(problem, {"dt":step})
    lr2 = RK2(problem, {"dt":step})
    lr3 = RK4(problem, {"dt":step})

    lr1.calculate()
    lr2.calculate()
    lr3.calculate()

    steps.append(step)
    error_r1.append(lr1.error(e))
    error_r2.append(lr2.error(e))
    error_r3.append(lr3.error(e))


matplotlib.pyplot.plot(steps, error_r1,'o-b',label="Euler")
matplotlib.pyplot.plot(steps, error_r2, 'd-r',label="Heun")
matplotlib.pyplot.plot(steps, error_r3, 'h-g',label="Runge-Kutta")

matplotlib.pyplot.yscale('log')
matplotlib.pyplot.ylabel('Error')

matplotlib.pyplot.xscale('log')
matplotlib.pyplot.xlabel('Step width')

matplotlib.pyplot.legend()
matplotlib.pyplot.show()