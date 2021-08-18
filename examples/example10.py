#
import sys

sys.path.append('../vodes')
#
 
from vodes.error.analysis import IntervalAnalysis as IA
from vodes.symbolic.expressions.rational import Rational

from pymbolic import var
from pymbolic.primitives import Power

##
# Example is from https://www.moodle.tum.de/pluginfile.php/2367064/mod_resource/content/3/NumPro_Vorlesung_Kapitel_1.pdf
##


# x² + 2px - q = 0
q = var('q')
p = var('p')

# 1.) x := sqrt(p²+q) - p
f1 = Power(p**2 + q,Rational(1,2)) - p
ia1 = IA(f1)

# 2.) x:= q / ( sqrt(p²+q)+p )
f2 = q / (Power(p**2 + q,Rational(1,2)) + p)
ia2 = IA(f2)

# p = 500, q = 1
context = {
    'p' : 500,
    'q' : 1
}

err1 = ia1.absolute(context=context, min_precision=11, max_precision=113)
err2 = ia2.absolute(context=context, min_precision=11, max_precision=113)

print(err1[0])
print(err2[0])

ia1.show()
ia2.show()