from numpy.lib.arraysetops import isin
from sympy.core import expr
from vodes.error.analysis import Analysis
from vodes.error.node import Node
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
        print("AFFINE : ")
        pprint(self._affine)

        self._approx = self._taylor_approximation(self._affine)
        print("APPROX : ")
        pprint(self._approx, use_unicode=True)

    # Responsible for injecting noise variables
    def _pre_walk(self, node, func, index=1):
        func(node, index)

        for i in range(0,len(node.args)):
            self._pre_walk(node.args[i], func, index * 10 + i)

    def _after_walk(self, node, func, index=1):
        for i in range(0,len(node.args)):
            self._after_walk(node.args[i], func, index * 10 + i)
        
        func(node, index)

    def _convert_to_tree(self, expr, tree=None):
        if tree is None:
            tree = Node(expr.func)
        
        for arg in expr.args:
            self._convert_to_tree(arg,tree=tree.insert(arg))

        return tree

    def _split(self, node, index):
        if len(node.args) <= 2:
            return node

        children = [
            node.args[0],
            node.args[1]
        ]

        for i in range (2, len(node.args)):
            children[0] = Node(node.func, [children[0],children[1]])
            children[1] = node.args[i]

        node.args = children  
        
    def _inject_noise(self, node, index):
        if len(node.args) == 0:
            return node

        _t = Node(node.func, node.args)
        noise = symbols(str(self._machine_epsilon) + str(index))

        node.func = Mul
        node.args = [
            _t,
           Node(Add, [
               Node(1),
               Node(noise)
           ]) 
        ]

        self._noises.append(noise)
        self._accurate[noise] = 0
        self._interval[noise] = SetExpr(Interval(-self._machine_epsilon,self._machine_epsilon))

        
    def _convert_to_binary(self, expr):
        tree = self._convert_to_tree(expr)
        self._pre_walk(tree, lambda n,i : self._split(n, i))
        return tree
  
    def _inject(self, expression):
        tree = self._convert_to_binary(expression)
        self._after_walk(tree, lambda n,i : self._inject_noise(n,i))
        return tree.get_expression()

    def _cleanup(self, expression):
        return self._interval_subs(expression,self._interval)

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
        #t_1 = reduce(lambda x,y: x + y, [diff(affine, noise).subs(self._accurate) * noise for noise in self._noises])
        t_1 = reduce(lambda x,y: x + y, [diff(affine, noise).subs(self._accurate) for noise in self._noises])
        approx_t_1 = self._machine_epsilon * Abs(t_1)

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
        #r_1 = reduce(lambda x,y: x + y,[diff(affine,e_i,e_j) * e_i * e_j for e_i in self._noises for e_j in self._noises])
        r_1 = reduce(lambda x,y: x + y,[diff(affine,e_i,e_j) for e_i in self._noises for e_j in self._noises])
        approx_r_1 = self._cleanup(self._machine_epsilon ** 2 * Abs(r_1))

        ##
        # Why we can not use subs on the whole expression, when usionvert=convertng interval arithmetic.
        # e,e0,e1 = symbols("e e0 e1", real=True, positive=True)
        # l = SetExpr(Interval(-e,e))
        # l * l = SetExpr(Interval(-e**2, e**2))
        # (e0 * e1).subs({e0:l,e1:l}) = SetExpr(Interval(0, e**2))
        ##

        return approx_t_1 + approx_r_1

    def _interval_add(self, sexpr, arg):
        smin = sexpr.set.start
        smax = sexpr.set.end

        amin = None
        amax = None 

        # ARG is SetExpr
        if isinstance(arg, SetExpr):
            amin = arg.set.start
            amax = arg.set.end
        else:
            amin = amax = arg

        return SetExpr(
            Interval(
                smin + amin,
                smax + amax
            )
        )

    def _interval_mul(self, sexpr, arg):
        smin = sexpr.set.start
        smax = sexpr.set.end

        amin = None
        amax = None 

        # ARG is SetExpr
        if isinstance(arg, SetExpr):
            amin = arg.set.start
            amax = arg.set.end
        else:
            amin = amax = arg
        
        return SetExpr(
            Interval(
                min(smin * amin, smin * amax, smax * amin, smax * amax),
                max(smin * amin, smin * amax, smax * amin, smax * amax)
            )
        )

    def _interval_calc(self, sexpr, sargs, expr):
        sres = sexpr[0]
        func = lambda res, arg : res.func(arg)

        if expr.is_Add:
            func = lambda res, arg : self._interval_add(res, arg)
        elif expr.is_Mul:
            func = lambda res, arg : self._interval_mul(res, arg)
        elif expr.is_Pow:
            func = lambda res, arg : res.__pow__(arg)
        ##
        # Inequality, boundary solution
        ##
        elif isinstance(expr, Abs):
            func = lambda res, arg : Abs(res.set.end)

        ##
        # Function using the set as only argument
        ##
        if len(expr.args) == 1:
            return func(sres, sres)

        for i in range(1,len(sexpr)):
            sres = func(sres, sexpr[i])

        missing = []
        for arg in sargs:
            try:
                sres = func(sres, arg)
            except TypeError:
                missing.append(arg)

        if len(missing) == 0:
            return sres
        else:
            missing.append(sres)
            return expr.func(*missing)

    def _interval_arithmetic(self, expr, args):
        res = None

        # Determine if interval contained
        sexpr = []
        sargs = []
        for c in args:
            if isinstance(c,SetExpr):
                sexpr.append(c)
            else:
                sargs.append(c)

        if len(sexpr) == 0:
            res = expr.func(*args) 
        else:
            res = self._interval_calc(sexpr, sargs, expr)
        
        return res
    
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
            child = expression.subs(subs)
            return child
        else:
            children = [self._interval_subs(x,subs) for x in expression.args]

        return self._interval_arithmetic(expression, children)
    
    # precision excludes! implicit bit
    def absolute(self, subs, precision=23):
        subs[self._machine_epsilon] = 2**(-precision)
        return self._interval_subs(self._approx, subs)