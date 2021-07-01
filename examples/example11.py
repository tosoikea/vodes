#
import sys
sys.path.append('../vodes')
#
 
from pymbolic import var
from vodes.error.analysis import IntervalAnalysis as IA

x = var("x")
y = var("y")
u = x * y

context ={
    "x" : 10,
    "y" : 20
}

ia = IA(problem=u)
err = ia.absolute(context=context)

