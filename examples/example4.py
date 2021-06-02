#
import sys
sys.path.append('../vodes')
#

from vodes.error import roundoff as err
from sympy import *

x = symbols('x')
y = symbols('y')
f = (x+y) * 5

e = err.Roundoff(f)
abs = e.absolute()

print(abs)