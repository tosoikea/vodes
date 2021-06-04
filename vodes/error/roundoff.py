from numpy.lib.arraysetops import isin
from vodes.error.analysis import Analysis
from abc import ABC, abstractmethod
from functools import reduce
# symbols, diff
from sympy import *
from sympy.sets.setexpr import SetExpr

class Roundoff(Analysis,ABC):
    def __init__(self, problem):
        self._problem = problem
        # contains all noise variables e_i
        self._noises = []
        # dictionary to allow substitution of all noise variables to 0 i.e. exact calculation
        self._accurate = {}
        self._interval = {}
        
        # with lowest precision (mantissa = 0) we can still display 1 (2^0 * 1) => e <= 1
        # it is impossible to exactly display 0 with an arbitrary amount of floating point precision => e > 0
        # => e \in (0,1]
        self._machine_epsilon = symbols("e", real=True, positive=True)

        self._affine = self._inject(self._problem)
        print("AFFINE : " + str(self._affine))
        self._approx = self._taylor_approximation(self._affine)
        print("APPROX : " + str(self._approx))

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
            self._interval[noise] = SetExpr(Interval(-self._machine_epsilon,self._machine_epsilon))
            # evaluate=False more in a debug sense, to allow for later validation of injection
            return Mul(expression.func(*[self._inject(x, index=index+str(ind)) for ind, x in enumerate(expression.args)],evaluate=False),(1 + noise),evaluate=False)

    def _cleanup(self, expression):
        return self._interval_subs(expression,{})

    def _taylor_approximation(self, affine):
        ##
        # Absolute Error := | o(f(x)) - f(x) |
        # Affine form with noise variables to model rounded operations
        #                <= | f(x,e) - f(x) |
        # Taylor relative to e at e = (0,...,0)
        # T_1 : first order term
        # R_1 : first order error term
        #                <= | f(x,0) + T_1(f,x,0) + R_1(f,x,p) - f(x) |
        # f(x,0) = f(x)
        #                 = | T_1(f,x,0) + R_1(f,x,p) |
        ##
        t_1 = Abs(reduce(lambda x,y: x + y, [diff(affine, noise) for noise in self._noises]).subs(self._accurate)) * self._machine_epsilon

        ##
        # Taylor’s Inequality :
        # M >= | f^(n+1)(x) | for all x \in (- \eps, + \eps) => | R_n(x) | <= M / (n+1)! * | x - a |^(n+1)
        #
        # R_1(x) = (x - a)^T / 2 * Hf(a)(x - a)
        #
        # d^2 f / de^2  = \sum_{i=0,j=0} (\partial^2 f) / (\partial e_i \partial e_j) (x,e)
        # https://mathinsight.org/taylors_theorem_multivariable_introduction
        #
        # Idea : Exact Remainder, Compute Intervall for noises, upper bound is M
        ##
        r_1 = Abs(reduce(lambda x,y: x + y,[diff(affine,e_i,e_j) for e_i in self._noises for e_j in self._noises]))
        
        return t_1 + r_1.subs(self._interval)
    
    def _interval_subs(self,expression,subs):
        ###
        # Unfortunately, sympy is displaying differentiating behavior for the following scenarios :
        # 1. 1 + SetExpr(Interval(-1.19209289550781E-7, 1.19209289550781E-7))
        #       = SetExpr(Interval(0.999999880790710, 1.00000011920929))
        # 2. (1 + SetExpr(Interval(-e, e))).subs(e,1.19209289550781E-7) 
        #       = SetExpr(ImageSet(Lambda(_d, _d + 1), Interval(-1.19209289550781E-7, 1.19209289550781E-7)))
        # 3. Add(1,SetExpr(Interval(-1.19209289550781E-7, 1.19209289550781E-7)))
        #       = 1 + SetExpr(Interval(-1.19209289550781E-7, 1.19209289550781E-7))
        #    ~ Tree traversal with leaf subs + e.func(*args) construction
        #
        # Note, that only the first version is correct and allows for future evaluation. ALL other cases lead to errors later on.
        # It seems only this version uses the operator functions defines within SetExpr while all others at some point default to the ImageSet. (lambda vs symbolic)
        # 
        # TODO : This bug has to be fixed on the side of sympy. -> Determine issue, fix bug, open PR
        # => Currently adds heavy overhead
        ###

        children = []
        if len(expression.args) == 0:
            return expression.subs(subs)
        else:
            children = [self._interval_subs(x,subs) for x in expression.args]

        sexpr = None
        sargs = []
        for c in children:
            if isinstance(c, SetExpr):
                sexpr = c
            else:
                sargs.append(c)
        
        if sexpr is None:
            res = expression.func(*children)
        else:
            if expression.is_Add:
                print("T1")
                res = sexpr.__add__(Add(*sargs))
            elif expression.is_Mul:
                print("T2")
                res = sexpr.__mul__(Mul(*sargs))
                print(res)
            elif expression.is_Pow:
                print("T3")
                res = sexpr.__pow__(*sargs)
            ##
            # Inequality, boundary solution
            ##
            elif isinstance(expression, Abs):
                if len(children) > 1:
                    raise Exception("Interval was not compacted")
                res = Abs(sexpr.set.end)
            else:
                res = expression.func(*children)

        return res
    
    # precision excludes! implicit bit
    def absolute(self, subs, precision=23):
        subs[self._machine_epsilon] = 2**(-precision)
        return self._interval_subs(self._approx, subs)