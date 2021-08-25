#
import sys

sys.path.append('../vodes')
#
 
from vodes.error.analysis import IntervalAnalysis as IA
from vodes.symbolic.expressions.nthroot import NthRoot

from pymbolic import var

##
# Example is from Introduction to Numerical Analysis from Arnold Neumaier (Cambridge) p.29
##

# Does not work because of https://github.com/sympy/sympy/issues/20097
# Wait for version 1.9 and replay test

x = var('x')

# 1.) p := sqrt(x+1) - sqrt(x)
f1 = NthRoot(x+1,n=2) - NthRoot(x,n=2)
ia1 = IA(f1)

# 2.) p := 1 / (sqrt(x+1) + sqrt(x))
f2 = 1 / (NthRoot(x+1,n=2) + NthRoot(x,n=2))
ia2 = IA(f2)

# x >> 1
context = {
    'x' : 420
}

err1 = ia1.absolute(context=context, min_precision=11, max_precision=53)
ia1.show()

err2 = ia2.absolute(context=context, min_precision=11, max_precision=53)
ia2.show()