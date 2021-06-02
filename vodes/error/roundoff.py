from vodes.error.analysis import Analysis
from abc import ABC, abstractmethod
from functools import reduce
# symbols, diff
from sympy import *

class Roundoff(Analysis,ABC):
    def __init__(self, problem):
        self._problem = problem
        # contains all noise variables e_i
        self._noises = []
        # dictionary to allow substitution of all noise variables to 0 i.e. exact calculation
        self._accurate = {}
        #
        self._machine_epsilon = symbols("e")
        
    # Responsible for injecting noise variables
    def _inject(self, expression, index="0"):
        # "Either it has empty args, in which case it is a leaf node in any expression tree"
        if len(expression.args) == 0:
            return expression
        # "or it has args, in which case, it is a branch node of any expression tree"
        else:
            noise = symbols(str(self._machine_epsilon) + index)
            self._noises.append(noise)
            self._accurate[noise] = 0
            return expression.func(*[self._inject(x, index=index+str(ind)) for ind, x in enumerate(expression.args)]) * (1 + noise)

    def _approximation(self, affine):
        # \eps * | first order taylor |
        first_order_approx = Abs(reduce(lambda x,y: x + y, [diff(affine, noise) for noise in self._noises]).subs(self._accurate)) * self._machine_epsilon
        #second_order_error = 0
        return first_order_approx
    
    def absolute(self):
        affine = self._inject(self._problem)
        print(affine)
        return self._approximation(affine)


