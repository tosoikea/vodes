#
import sys
sys.path.append('../vodes')
#

from vodes.ode.problem import Problem
from vodes.ode.runge import RK1, RK2
from vodes.symbolic.expressions.rational import Rational
from vodes.error.analysis import IntervalAnalysis as IA, TaylorAnalysis as TA
from vodes.error.utils import AnalysisSolution, PseudoExactODESolution, show
from vodes.symbolic.mapper.extended_substitution_mapper import substitute

from pymbolic.primitives import Variable

MIN_PREC = 24
MAX_PREC = 53

MIN_EXP = 8
MAX_EXP = 8

iv = Rational(1,2)
dt1 = Rational(1,3)
dt2 = Rational(1,5)

# define function
t = Variable('t')
y = Variable('y')
f = -4 * y * t

problem = Problem(f, iv, (0,3))

lr1 = RK1(problem, {"dt":dt1})
lr2 = RK2(problem, {"dt":dt1})

##
p1_expr = None
p2_expr = None

for i in range(0,len(lr1.exprs)):
    (_,expr) = lr1.exprs[i]
    if p1_expr is None:
        p1_expr = expr
    else:
        p1_expr = substitute(
            expr, {
                f'y_{i-1}' : p1_expr
            }
        )

for i in range(0,len(lr2.exprs)):
    (_,expr) = lr2.exprs[i]
    if p2_expr is None:
        p2_expr = expr
    else:
        p2_expr = substitute(
            expr, {
                f'y_{i-1}' : p2_expr
            }
        )

ta1 = TA(p1_expr)
ta2 = TA(p2_expr)

context = {}
errt1 = ta1.absolute(context=context, min_precision=MIN_PREC, max_precision=MAX_PREC, min_exponent=MIN_EXP,max_exponent=MAX_EXP)
errt2 = ta2.absolute(context=context, min_precision=MIN_PREC, max_precision=MAX_PREC, min_exponent=MIN_EXP,max_exponent=MAX_EXP)

show(
    solutions=[
        AnalysisSolution(bexprs=errt1,name="TA(RK1)"),
        AnalysisSolution(bexprs=errt2,name="TA(RK2)"),
        PseudoExactODESolution(solver=lr1,name="RK1"),
        PseudoExactODESolution(solver=lr2,name="RK2")
    ],
    min_prec=MIN_PREC,
    max_prec=MAX_PREC
)