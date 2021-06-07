#
import sys
sys.path.append('../vodes')
#

from vodes.error import roundoff as err
from sympy import *

a = symbols('a')
h = symbols('h', positive=True)
c = 3.0

##
# DO NOT USE x y!
##

f = c * a
# euler
step = a + h * f

#step = step.subs(h,0.5)
print("EXAMPLE 4")
print("---------")
pprint(step)

e = err.Roundoff(step)
abs = e.absolute({a : 0.000001, h : 0.025})
#abs = e.absolute({})

print(abs)