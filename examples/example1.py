#
import sys
sys.path.append('../vodes')
#

import vodes
import vodes.runge
import math

iv = 1
c = 1.0
f = lambda t,y : c * y
e = lambda t : math.exp(c * t) * iv

lr1 = vodes.runge.RK1(f, iv, (0,1), {"dt":0.5})
lr1.calculate()
lr1.show(exact=e)

lr2 = vodes.runge.RK2(f, iv, (0,1), {"dt":0.001})
lr2.calculate()
lr2.show(exact=e)