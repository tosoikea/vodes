#
import sys
sys.path.append('../vodes')
#

from vodes.ode.problem import Problem
from vodes.ode.runge import RK1
from sympy import *

iv = 1/2

# define function
t = symbols('t')
y = symbols('y')
f = -4 * y * t

problem = Problem(f, iv, (0,5))
lr1 = RK1(problem, {"dt":0.1})
lr1.calculate()
lr1.show()

lr2 = RK1(problem, {"dt":0.001})
lr2.calculate()
lr2.show()