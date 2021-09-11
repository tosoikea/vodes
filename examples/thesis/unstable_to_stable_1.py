#
import sys

sys.path.append('../vodes')
#
 
from vodes.error.roundoff_analysis import IntervalAnalysis as IA, TaylorAnalysis as TA
from vodes.symbolic.expressions.nthroot import NthRoot
from vodes.symbolic.expressions.primitives import Subtraction

from pymbolic import var
from pymbolic.primitives import Quotient
from vodes.error.utils import PseudoExactSolution, show, export, AnalysisSolution
from mpmath import mpf, sqrt

MIN_PREC = 11
MAX_PREC = 53

MIN_EXP = 8
MAX_EXP = 8

##
# Example is from Introduction to Numerical Analysis from Arnold Neumaier (Cambridge) p.29
##

x = var('x')
xv = Quotient(1277,10)

# 1.) p := sqrt(x+1) - sqrt(x)
f1 = Subtraction((NthRoot(x+1,n=2),NthRoot(x,n=2)))
ia1 = IA(f1)
ta1 = TA(f1)
ex1 = lambda: sqrt(mpf('127.7') + mpf(1)) - sqrt(mpf('127.7'))

# 2.) p := 1 / (sqrt(x+1) + sqrt(x))
f2 = 1 / (NthRoot(x+1,n=2) + NthRoot(x,n=2))
ia2 = IA(f2)
ta2 = TA(f2)
ex2 = lambda: mpf(1) / (sqrt(mpf('127.7') + mpf(1)) + sqrt(mpf('127.7')))

# x >> 1
context = {
    'x' : xv
}

erri1 = ia1.absolute(context=context, min_precision=MIN_PREC, max_precision=MAX_PREC, min_exponent=MIN_EXP, max_exponent=MAX_EXP)
erri2 = ia2.absolute(context=context, min_precision=MIN_PREC, max_precision=MAX_PREC, min_exponent=MIN_EXP, max_exponent=MAX_EXP)
errt1 = ta1.absolute(context=context, min_precision=MIN_PREC, max_precision=MAX_PREC, min_exponent=MIN_EXP, max_exponent=MAX_EXP)
errt2 = ta2.absolute(context=context, min_precision=MIN_PREC, max_precision=MAX_PREC, min_exponent=MIN_EXP, max_exponent=MAX_EXP)

export(
    file="unstable_to_stable_1.csv",
    solutions=[
        AnalysisSolution(bexprs=erri1,name="IA_SUB"),
        AnalysisSolution(bexprs=erri2,name="IA_DIV"),
        AnalysisSolution(bexprs=errt1,name="TA_SUB"),
        AnalysisSolution(bexprs=errt2,name="TA_DIV"),
        PseudoExactSolution(func=ex1,name="Exact_SUB"),
        PseudoExactSolution(func=ex2,name="Exact_DIV")
    ],
    min_prec=MIN_PREC,
    max_prec=MAX_PREC
)