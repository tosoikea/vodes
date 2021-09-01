#
import sys

sys.path.append('../vodes')
#
 
from vodes.error.analysis import IntervalAnalysis as IA, TaylorAnalysis as TA
from vodes.symbolic.expressions.nthroot import NthRoot
from vodes.symbolic.expressions.primitives import Subtraction

from pymbolic import var
from vodes.error.utils import PseudoExactSolution, show, AnalysisSolution
from mpmath import mpf, sqrt

MIN_PREC = 11
MAX_PREC = 53

##
# Example is from Introduction to Numerical Analysis from Arnold Neumaier (Cambridge) p.29
##

# Does not work with TA symbolic because of https://github.com/sympy/sympy/issues/20097

x = var('x')
xv = 420

# 1.) p := sqrt(x+1) - sqrt(x)
f1 = Subtraction((NthRoot(x+1,n=2),NthRoot(x,n=2)))
#ia1 = IA(f1)
#ta1 = TA(f1)
ex1 = lambda: sqrt(mpf(xv) + mpf(1)) - sqrt(mpf(xv))

# 2.) p := 1 / (sqrt(x+1) + sqrt(x))
f2 = 1 / (NthRoot(x+1,n=2) + NthRoot(x,n=2))
#ia2 = IA(f2)
#ta2 = TA(f2)
ex2 = lambda: mpf(1) / (sqrt(mpf(xv) + mpf(1)) + sqrt(mpf(xv)))

# x >> 1
context = {
    'x' : xv
}

#erri1 = ia1.absolute(context=context, min_precision=MIN_PREC, max_precision=MAX_PREC)
#errt1 = ta1.absolute(context=context, min_precision=MIN_PREC, max_precision=MAX_PREC)

#erri2 = ia2.absolute(context=context, min_precision=MIN_PREC, max_precision=MAX_PREC)
#errt2 = ta2.absolute(context=context, min_precision=MIN_PREC, max_precision=MAX_PREC)

show(
    solutions=[
        #AnalysisSolution(bexprs=erri1,name="IA (SUB)"),
        #AnalysisSolution(bexprs=errt1,name="TA (SUB)"),
        PseudoExactSolution(func=ex1,name="Exact (SUB)"),
        #AnalysisSolution(bexprs=erri2,name="IA (DIV)"),
        #AnalysisSolution(bexprs=errt2,name="TA (DIV)"),
        #PseudoExactSolution(func=ex2,name="Exact (DIV)")
    ],
    min_prec=MIN_PREC,
    max_prec=MAX_PREC
)