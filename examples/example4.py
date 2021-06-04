#
import sys
sys.path.append('../vodes')
#

from vodes.error import roundoff as err
from sympy import *

a = symbols('a')

##
# DO NOT USE x y!
##


f = (a+1) * 1/a
print(f)
e = err.Roundoff(f)
abs = e.absolute({a : 0.01})
#abs = e.absolute({})

print(abs)