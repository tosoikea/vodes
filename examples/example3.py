#
import sys
sys.path.append('../vodes')
#

from vodes import runge
from sympy import symbols
from sympy.functions import exp

iv = 1
c = 1.0

# define function
t = symbols('t')
y = symbols('y')
f = t * y

# exact solution
e = exp(c * t) * iv

lr1 = runge.RK1(f, iv, (0,1), {"dt":0.5}, symbols=[y, t])
lr1.calculate()
lr1.show(exact=e)

lr2 = runge.RK2(f, iv, (0,1), {"dt":0.001}, symbols=[y, t])
lr2.calculate()
lr2.show(exact=e)