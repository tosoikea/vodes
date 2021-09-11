#
import sys

sys.path.append('../vodes')
#
 
from vodes.symbolic.expressions.primitives import Subtraction
from vodes.symbolic.expressions.nthroot import NthRoot
from vodes.error.roundoff_analysis import IntervalAnalysis as IA, TaylorAnalysis as TA

from vodes.error.utils import PseudoExactSolution, export, show, AnalysisSolution
from vodes.symbolic.mapper.extended_evaluation_mapper import evaluate

from mpmath import mpf, sqrt
from pymbolic import var
from pymbolic.primitives import Quotient
from vodes.symbolic.expressions.interval import Interval

MIN_PREC = 11
MAX_PREC = 53

MIN_EXP = 8
MAX_EXP = 8

##
# Example is from https://fpbench.org/benchmarks.html
##

x1 = var('x1')
x2 = var('x2')
sym_xv1 = Interval(1,5)
sym_xv2 = Interval(1,5)

xv1 = evaluate(sym_xv1)
xv2 = evaluate(sym_xv2)

a = Subtraction(((x1 * x1) + x2,11))
b = Subtraction((x1 + (x2 * x2),7))

p = ((a*a) + (b*b))

ta2 = TA(p)

# x ~ 0
context = {
    'x1' : sym_xv1,
    'x2' : sym_xv2
}

errt2 = ta2.absolute(context=context, min_precision=MIN_PREC, max_precision=MAX_PREC, min_exponent=MIN_EXP,max_exponent=MAX_EXP)

export(
    file='himmilbeau_ranged_analysis.csv',
    solutions=[
        AnalysisSolution(bexprs=errt2,name="TA")
    ],
    min_prec=MIN_PREC,
    max_prec=MAX_PREC
)