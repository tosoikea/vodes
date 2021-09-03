#
import sys
sys.path.append('../vodes')
#

from vodes.ode.problem import Problem
from vodes.ode.runge import RK1
from vodes.symbolic.expressions.rational import Rational

from pymbolic.primitives import Variable

PREC = 53

iv = Rational(1,2)
dt1 = Rational(1,5)
dt2 = Rational(1,10)

# define function
t = Variable('t')
y = Variable('y')
f = -4 * y * t

problem = Problem(f, iv, (0,3))
lr1 = RK1(problem, {"dt":dt1})
lr1.calculate(precision=PREC)
lr1.show()

lr2 = RK1(problem, {"dt":dt2})
lr2.calculate(precision=PREC)
lr2.show()