#
import sys

sys.path.append('../vodes')
#
 
from vodes.symbolic.expressions.rational import Rational
from vodes.error.analysis import IntervalAnalysis as IA

from pymbolic import var
from pymbolic.primitives import Power

##
# Example is from https://www.math.utk.edu/~ccollins/M577/Handouts/cond_stab.pdf
##


x = var('x')

# 1.) x := sqrt(1+x) - 1
f1 = Power(1+x,Rational(1,2)) - 1
print(f1)
ia1 = IA(f1)

# 2.) x:= x / (sqrt(1+x) + 1)
f2 = x / (Power(1+x,Rational(1,2)) + 1)
ia2 = IA(f2)

# x = 0
context = {
    'x' : 0
}

err1 = ia1.absolute(context=context, min_precision=11, max_precision=113)
print(err1[0])
ia1.show()

err2 = ia2.absolute(context=context, min_precision=11, max_precision=113)
print(err2[0])
ia2.show()