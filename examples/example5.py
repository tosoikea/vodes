#
import sys
sys.path.append('../vodes')
#

from vodes.ode.problem import Problem
from vodes.ode.runge import RK1, RK2, RK4
from sympy import *

# define function
t = symbols('t')
y = MatrixSymbol('y',2,1)

iv = Matrix([1/2, -1/3])
f = Matrix([[0,-2],[4,0]]) * y

problem = Problem(f, iv, (0,50))
lr1 = RK1(problem, {"dt":0.01})
lr1.calculate()
lr1.show()

lr2 = RK2(problem, {"dt":0.01})
lr2.calculate()
lr2.show()

lr3 = RK4(problem, {"dt":0.01})
lr3.calculate()
lr3.show()