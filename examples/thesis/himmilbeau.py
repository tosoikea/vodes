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

MIN_PREC = 11
MAX_PREC = 53

MIN_EXP = 8
MAX_EXP = 8

##
# Example is from https://fpbench.org/benchmarks.html
##

x1 = var('x1')
x2 = var('x2')
sym_xv1 = Quotient(1,10)
sym_xv2 = Quotient(3,10)

xv1 = evaluate(sym_xv1)
xv2 = evaluate(sym_xv2)

a = (((x1 * x1) + x2) - 11)
b = ((x1 + (x2 * x2)) - 7)

p = ((a*a) + (b*b))

ia1 = IA(p)
ta1 = TA(p,{"n":1})
ta2 = TA(p,{"n":2})

ax1 = lambda: mpf(xv1) * mpf(xv2) + mpf(xv2) - mpf(11)
bx1 = lambda: mpf(xv1) + mpf(xv2) * mpf(xv2) - mpf(7)
ex1 = lambda: ax1() * ax1() + bx1() * bx1()

# x ~ 0
context = {
    'x1' : sym_xv1,
    'x2' : sym_xv2
}

erri1 = ia1.absolute(context=context, min_precision=MIN_PREC, max_precision=MAX_PREC, min_exponent=MIN_EXP,max_exponent=MAX_EXP)
errt1 = ta1.absolute(context=context, min_precision=MIN_PREC, max_precision=MAX_PREC, min_exponent=MIN_EXP,max_exponent=MAX_EXP)
errt2 = ta2.absolute(context=context, min_precision=MIN_PREC, max_precision=MAX_PREC, min_exponent=MIN_EXP,max_exponent=MAX_EXP)

show(
    solutions=[
        AnalysisSolution(bexprs=erri1,name="IA"),
        AnalysisSolution(bexprs=errt1,name="TA(n=1)"),
        AnalysisSolution(bexprs=errt2,name="TA(n=2)"),
        PseudoExactSolution(func=ex1),
    ]
)