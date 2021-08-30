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
# Example is from https://www.math.utk.edu/~ccollins/M577/Handouts/cond_stab.pdf
##


x = var('x')
sym_xv = Quotient(1,10*6)
xv = evaluate(sym_xv)

# 1.) x := sqrt(1+x) - 1
f1 = NthRoot(1+x,2) - 1
ia1 = IA(f1)
ta1 = TA(f1)
ex1 = lambda: sqrt(mpf(1) + mpf(xv)) - mpf(1)

# 2.) x:= x / (sqrt(1+x) + 1)
f2 = x / (NthRoot(1+x,2) + 1)
ia2 = IA(f2)
ta2 = TA(f2)
ex2 = lambda: mpf(xv) / (sqrt(mpf(1) + mpf(xv)) + mpf(1))

# x ~ 0
context = {
    'x' : sym_xv
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

erri2 = ia2.absolute(context=context, min_precision=11, max_precision=53)
errt2 = ta2.absolute(context=context, min_precision=11, max_precision=53)
show(
    solutions=[
        AnalysisSolution(bexprs=erri2,name="IA"),
        AnalysisSolution(bexprs=errt2,name="TA"),
        PseudoExactSolution(func=ex2),
    ]
)
