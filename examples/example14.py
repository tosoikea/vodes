#
import sys

sys.path.append('../vodes')
#
 
from vodes.symbolic.expressions.primitives import Subtraction
from vodes.error.analysis import IntervalAnalysis as IA, TaylorAnalysis as TA

from vodes.error.utils import PseudoExactSolution, show, AnalysisSolution
from vodes.symbolic.mapper.extended_evaluation_mapper import evaluate

from mpmath import mpf, sqrt
from pymbolic import var
from pymbolic.primitives import Quotient

##
# Example is from thesis
##

x = var('x')
xv = 800

# 1.) x := sqrt(1+x) - 1
f1 = (Subtraction((x,x)))**20
ia1 = IA(f1)
ta1 = TA(f1)
ex1 = lambda: sqrt(mpf(xv) + mpf(xv)) - mpf(1)

# x ~ 0
context = {
    'x' : xv
}

erri1 = ia1.absolute(context=context, min_precision=11, max_precision=53)
errt1 = ta1.absolute(context=context, min_precision=11, max_precision=53)
show(
    solutions=[
        AnalysisSolution(bexprs=erri1,name="IA"),
        AnalysisSolution(bexprs=errt1,name="TA"),
        PseudoExactSolution(func=ex1),
    ]
)