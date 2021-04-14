#
import sys
sys.path.append('../vodes')
#

import vodes
import vodes.runge
import math
import matplotlib.pyplot

iv = 1
c = 3.0
f = lambda t,y : c * y
e = lambda t : math.exp(c * t) * iv

# visualize
steps = []
error_r1 = []
error_r2 = []

for i in range(0,8):
    step = 10 ** -i
    print(step)
    lr1 = vodes.runge.RK1(f, iv, (0,1), {"dt":step})
    lr2 = vodes.runge.RK2(f, iv, (0,1), {"dt":step})

    lr1.calculate()
    lr2.calculate()

    steps.append(step)
    error_r1.append(lr1.error(e))
    error_r2.append(lr2.error(e))


matplotlib.pyplot.plot(steps, error_r1,'o-b',label="Euler")
matplotlib.pyplot.plot(steps, error_r2, 'd-r',label="Heun")

matplotlib.pyplot.yscale('log')
matplotlib.pyplot.ylabel('Error')

matplotlib.pyplot.xscale('log')
matplotlib.pyplot.xlabel('Step width')

matplotlib.pyplot.legend()
matplotlib.pyplot.show()