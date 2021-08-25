#
import sys

sys.path.append('../vodes')
#
 
from vodes.symbolic.expressions.nthroot import NthRoot
from vodes.error.analysis import IntervalAnalysis as IA

from pymbolic import var
from pymbolic.primitives import Quotient

##
# Example is from https://www.math.utk.edu/~ccollins/M577/Handouts/cond_stab.pdf
##


x = var('x')

# 1.) x := sqrt(1+x) - 1
f1 = NthRoot(1+x,2) - 1
print(f1)
ia1 = IA(f1)

# 2.) x:= x / (sqrt(1+x) + 1)
f2 = x / (NthRoot(1+x,2) + 1)
ia2 = IA(f2)

# x ~ 0
context = {
    'x' : Quotient(1,10*6)
}

err1 = ia1.absolute(context=context, min_precision=11, max_precision=113)
ia1.show()

err2 = ia2.absolute(context=context, min_precision=11, max_precision=113)
ia2.show()