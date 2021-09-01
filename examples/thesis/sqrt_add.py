#
import sys

sys.path.append('../vodes')
#
 
from vodes.symbolic.expressions.nthroot import NthRoot
from vodes.error.analysis import IntervalAnalysis as IA, TaylorAnalysis as TA

from vodes.error.utils import PseudoExactSolution, show, AnalysisSolution
from vodes.symbolic.mapper.extended_evaluation_mapper import evaluate

from mpmath import mpf, sqrt
from pymbolic import var
from pymbolic.primitives import Quotient

##
# Example is from https://fpbench.org/benchmarks.html
##

x = var('x')
xv = 10

# p := 1 / (sqrt(x+1) + sqrt(x))
f1 = 1 / (NthRoot(x+1,2) + NthRoot(x,2))
ia1 = IA(f1)
ta1 = TA(f1)
ex1 = lambda: mpf(1) / (sqrt(mpf(xv) + mpf(1)) + sqrt(mpf(xv)))

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