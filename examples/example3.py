#
import sys
sys.path.append('../vodes')
#

from vodes.ode.runge import RK1, RK2
from vodes.ode.problem import Problem
from sympy import symbols
from sympy.functions import exp

iv = 1
c = 1.0

# define function
t = symbols('t')
y = symbols('y')
f = t * y
problem = Problem(f, iv, (0,1))

# exact solution
e = exp(c * t) * iv

lr1 = RK1(problem, {"dt":0.5})
lr1.calculate()
lr1.show(exact=e)

lr2 = RK2(problem, {"dt":0.001})
lr2.calculate()
lr2.show(exact=e)