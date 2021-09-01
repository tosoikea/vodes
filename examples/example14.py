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

MIN_PREC = 11
MAX_PREC = 53

##
# Example is from thesis
##

x = var('x')
xv = 800

# 1.) x := (x-x)^5
f1 = (Subtraction((x,x)))**5
ia1 = IA(f1)
ta1 = TA(f1)
ex1 = lambda: (mpf(xv) - mpf(xv))**5

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
        PseudoExactSolution(func=ex1),
    ],
    min_prec=MIN_PREC,
    max_prec=MAX_PREC
)