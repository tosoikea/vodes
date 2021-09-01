#
import sys

sys.path.append('../vodes')
#
 
from vodes.error.analysis import IntervalAnalysis as IA, TaylorAnalysis as TA
from vodes.symbolic.expressions.nthroot import NthRoot

from pymbolic import var
from vodes.error.utils import show, AnalysisSolution

MIN_PREC = 11
MAX_PREC = 53

##
# Example is from https://www.moodle.tum.de/pluginfile.php/2367064/mod_resource/content/3/NumPro_Vorlesung_Kapitel_1.pdf
##

# x² + 2px - q = 0
q = var('q')
p = var('p')

# 1.) x := sqrt(p²+q) - p
f1 = NthRoot(p**2 + q,n=2) - p
ia1 = IA(f1)
ta1 = TA(f1)

# 2.) x:= q / ( sqrt(p²+q)+p )
f2 = q / (NthRoot(p**2 + q,n=2) + p)
ia2 = IA(f2)
ta2 = TA(f2)

# p = 500, q = 1
context = {
    'p' : 500,
    'q' : 1
}

erri1 = ia1.absolute(context=context, min_precision=MIN_PREC, max_precision=MAX_PREC)
errt1 = ta1.absolute(context=context, min_precision=MIN_PREC, max_precision=MAX_PREC)

erri2 = ia2.absolute(context=context, min_precision=MIN_PREC, max_precision=MAX_PREC)
errt2 = ta2.absolute(context=context, min_precision=MIN_PREC, max_precision=MAX_PREC)

show(
    solutions=[
        AnalysisSolution(bexprs=erri1,name="IA (SUB)"),
        AnalysisSolution(bexprs=errt1,name="TA (SUB)"),
        AnalysisSolution(bexprs=erri2,name="IA (DIV)"),
        AnalysisSolution(bexprs=errt2,name="TA (DIV)")
    ]
)
