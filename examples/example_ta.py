#
import sys

sys.path.append('../vodes')
#
 
from vodes.error.analysis import IntervalAnalysis as IA, TaylorAnalysis as TA
from vodes.symbolic.expressions.nthroot import NthRoot
from vodes.symbolic.expressions.primitives import Subtraction
from vodes.symbolic.expressions.rational import Rational

from pymbolic import var
from vodes.error.utils import PseudoExactSolution, show, AnalysisSolution
from mpmath import mpf, sqrt

MIN_PREC = 11
MAX_PREC = 32

MIN_EXP = 8
MAX_EXP = 8

##
# Example is from Introduction to Numerical Analysis from Arnold Neumaier (Cambridge) p.29
##

x = var('x')
y = var('y')

xv = 420
yv = Rational(1,10)

# 1.) p := sqrt(x+1) - sqrt(x)
f = NthRoot(x+1,2) - NthRoot(x,2)
ta1 = TA(f,opt={"n":1})
ta2 = TA(f,opt={"n":2})
ta3 = TA(f,opt={"n":3})

ex = lambda: sqrt(mpf(xv) - mpf(0.1)) + sqrt(mpf(0.1))**2

# x >> 1
context = {
    'x' : xv,
    'y' : yv
}

err1 = ta1.absolute(context=context, min_precision=MIN_PREC, max_precision=MAX_PREC, min_exponent=MIN_EXP, max_exponent=MAX_EXP)
err2 = ta2.absolute(context=context, min_precision=MIN_PREC, max_precision=MAX_PREC, min_exponent=MIN_EXP, max_exponent=MAX_EXP)
err3 = ta3.absolute(context=context, min_precision=MIN_PREC, max_precision=MAX_PREC, min_exponent=MIN_EXP, max_exponent=MAX_EXP)

show(
    solutions=[
        AnalysisSolution(bexprs=err1,name="TA(n=1)"),
        AnalysisSolution(bexprs=err2,name="TA(n=2)"),
        AnalysisSolution(bexprs=err2,name="TA(n=3)"),
        PseudoExactSolution(func=ex)
    ],
    min_prec=MIN_PREC,
    max_prec=MAX_PREC
)