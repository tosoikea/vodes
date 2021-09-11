#
import sys

sys.path.append('../vodes')
#
 
from vodes.error.roundoff_analysis import IntervalAnalysis as IA, TaylorAnalysis as TA
from vodes.symbolic.expressions.nthroot import NthRoot

from pymbolic import var
from pymbolic.primitives import Quotient
from vodes.error.utils import PseudoExactSolution, show, export, AnalysisSolution
from mpmath import mpf, sqrt

MIN_PREC = 11
MAX_PREC = 53

MIN_EXP = 8
MAX_EXP = 8

##
# Example is from https://www.moodle.tum.de/pluginfile.php/2367064/mod_resource/content/3/NumPro_Vorlesung_Kapitel_1.pdf
##

# x² + 2px - q = 0
q = var('q')
qv = Quotient(99,100)
p = var('p')
pv = 121 + Quotient(3,10)

# 1.) x := sqrt(p²+q) - p
f1 = NthRoot(p**2 + q,n=2) - p
ia1 = IA(f1)
ta1 = TA(f1)
ex1 = lambda: sqrt(mpf('121.3')**2 + mpf('0.99')) - mpf('121.3')

# 2.) x:= q / ( sqrt(p²+q)+p )
f2 = q / (NthRoot(p**2 + q,n=2) + p)
ia2 = IA(f2)
ta2 = TA(f2)
ex2 = lambda: mpf('0.99') / (sqrt(mpf('121.3')**2 + mpf('0.99')) + mpf('121.3'))

# p=121.3,q=0.99
context = {
    'p' : pv,
    'q' : qv
}

erri1 = ia1.absolute(context=context, min_precision=MIN_PREC, max_precision=MAX_PREC, min_exponent=MIN_EXP, max_exponent=MAX_EXP)
erri2 = ia2.absolute(context=context, min_precision=MIN_PREC, max_precision=MAX_PREC, min_exponent=MIN_EXP, max_exponent=MAX_EXP)
errt1 = ta1.absolute(context=context, min_precision=MIN_PREC, max_precision=MAX_PREC, min_exponent=MIN_EXP, max_exponent=MAX_EXP)
errt2 = ta2.absolute(context=context, min_precision=MIN_PREC, max_precision=MAX_PREC, min_exponent=MIN_EXP, max_exponent=MAX_EXP)

export(
    file="unstable_to_stable_2.csv",
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
