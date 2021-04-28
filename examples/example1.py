#
import sys
sys.path.append('../vodes')
#

from vodes import runge
from sympy import *

iv = frac(Rational(1,2))
c = 1.0

# define function
t = symbols('t')
y = symbols('y')
f = -4 * y

lr1 = runge.RK1(f, iv, (0,5), {"dt":0.1}, symbols=[y, t])
lr1.calculate()
lr1.show()

lr2 = runge.RK2(f, iv, (0,5), {"dt":0.001}, symbols=[y, t])
lr2.calculate()
lr2.show()