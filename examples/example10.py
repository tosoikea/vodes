#
import sys

sys.path.append('../vodes')
#
 
from pymbolic import var
from vodes.symbolic.interval import ExactIntersectionEvaluator as EIE, Interval
from vodes.symbolic.symbols import BoundedVariable, MachineError
from vodes.error.analysis import IntervalAnalysis as IA

##
# Example is from https://www.moodle.tum.de/pluginfile.php/2367064/mod_resource/content/3/NumPro_Vorlesung_Kapitel_1.pdf
##


# x² + 2px - q = 0
q = var('q')
p = var('p')

# 1.) x := sqrt(p²+q) - p
f1 = (p**2 + q)**(0.5) - p
ia1 = IA(f1)

# 2.) x:= q / ( sqrt(p²+q)+p )
f2 = q / ((p**2 + q)**(0.5) + p)
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