#
import sys

sys.path.append('../vodes')
#
 
from vodes.symbolic.expressions.primitives import Subtraction
from vodes.error.analysis import IntervalAnalysis as IA, TaylorAnalysis as TA

from vodes.error.utils import PseudoExactIntervalSolution, show, AnalysisSolution
from vodes.symbolic.mapper.extended_evaluation_mapper import evaluate
from vodes.symbolic.expressions.interval import Interval

from mpmath import iv
from pymbolic import var
from pymbolic.primitives import Quotient

MIN_PREC = 11
MAX_PREC = 16

##
# Simple Interval Input Example
##

x = var('x')
xv = Interval(3,4)

# f := x**2 - x

f1 = Subtraction((x**2,x))
ia1 = IA(f1)
ta1 = TA(f1)
ex1 = lambda: iv.mpf([3,4])**2 - iv.mpf([3,4])

# x ~ 0
context = {
    'x' : xv
}

erri1 = ia1.absolute(context=context, min_precision=MIN_PREC, max_precision=MAX_PREC)
errt1 = ta1.absolute(context=context, min_precision=MIN_PREC, max_precision=MAX_PREC)

show(
    solutions=[
        AnalysisSolution(bexprs=erri1,name="IA"),
        AnalysisSolution(bexprs=errt1,name="TA"),
        PseudoExactIntervalSolution(func=ex1),
    ],
    min_prec=MIN_PREC,
    max_prec=MAX_PREC
)