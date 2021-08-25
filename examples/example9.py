#
import sys

sys.path.append('../vodes')
#
 
from vodes.error.analysis import IntervalAnalysis as IA
from vodes.symbolic.expressions.nthroot import NthRoot
from vodes.symbolic.expressions.trigonometric import sin,cos

from pymbolic import var

##
# Test trigonometric functions
##

x = var('x')

f = sin(x + 1) * cos(x)
ia = IA(f)

context = {
    'x' : 250
}

err1 = ia.absolute(context=context, min_precision=11, max_precision=53)
ia.show()