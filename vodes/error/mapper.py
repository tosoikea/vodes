from functools import reduce
import logging

from vodes.symbolic.mapper.interop import ExactPymbolicToSympyMapper, ExactSympyToPymbolicMapper
from vodes.symbolic.expressions.nthroot import NthRoot
from vodes.symbolic.expressions.primitives import Subtraction
from vodes.symbolic.expressions.trigonometric import sin, cos
from vodes.symbolic.expressions.bounded import MachineError, SmallestMachineNumber
from vodes.symbolic.expressions.interval import Interval

from pymbolic.mapper import RecursiveMapper
from pymbolic.primitives import Power, Variable

class ErrorMapper(RecursiveMapper):
    def __init__(self,context:dict,min_precision:int,max_precision:int,min_exponent:int,max_exponent:int):
        if self is ErrorMapper:
            raise TypeError("ErrorAnalysis may not be directly instantiated")
        
        ## Substitution of free variables
        self.context = context

        ## Analysis information
        self.min_precision = min_precision
        self.max_precision = max_precision
        self.min_exponent = min_exponent
        self.max_exponent = max_exponent

    @property
    def sigm(self):
        return SmallestMachineNumber(min_exponent=self.min_exponent,max_exponent=self.max_exponent)

    @property
    def eps(self):
        return MachineError(min_precision=self.min_precision,max_precision=self.max_precision)
        
    def map_variable(self, expr):
        from pymbolic.mapper.evaluator import UnknownVariableError

        try:
            subs = self.context[expr.name]
            
            # Obv. I can still create substitutions like a = b + 1, b = a + 1.
            # However, this is the users fault :)
            if isinstance(subs,Variable):
                raise ValueError("Invalid mapping of variable to variable.")

            return self.rec(subs)
        except KeyError:
            raise UnknownVariableError(expr.name)

    def map_rational(self, expr):
        from pymbolic.primitives import Quotient
        return self.map_quotient(Quotient(expr.numerator, expr.denominator))

    def map_call(self, expr):
        return self.rec(expr.function)(*[self.rec(par) for par in expr.parameters])
        

class IntervalMapper(ErrorMapper):
    def __init__(self,context:dict,min_precision:int,max_precision:int,min_exponent:int,max_exponent:int):
        super().__init__(context=context,min_precision=min_precision,max_precision=max_precision,min_exponent=min_exponent,max_exponent=max_exponent)

    def _error(self, expr, sigma:int, epsilon:int):
        # x(1+e) + d
        res = expr * Interval(1 - self.eps, 1 + self.eps)
        
        # TODO : Underflow
        #+ sigma * Interval(self.sigm.bound.start,self.sigm.bound.end)
        return res

    def map_constant(self, expr):
        if isinstance(expr,int):
            return expr
        
        return self._error(expr,0,1)

    def map_rational(self, expr):
        return self._error(expr,0,1)

    def map_sum(self, expr):
        from vodes.symbolic.mapper.extended_evaluation_mapper import evaluate
        assert len(expr.children) == 2

        l = self.rec(expr.children[0])
        r = self.rec(expr.children[1])

        # TODO : Define different strategies to detect subnormalities, epsilon differentiations etc.
        sigma = 1
        epsilon = 1

        if isinstance(sigma,int) and isinstance(epsilon,int):
            sigma = 0
            epsilon = 0

        return self._error(l+r,sigma=sigma,epsilon=epsilon)

    def map_sub(self, expr):
        from vodes.symbolic.mapper.extended_evaluation_mapper import evaluate
        assert len(expr.children) == 2

        l = self.rec(expr.children[0])
        r = self.rec(expr.children[1])

        # TODO : Define different strategies to detect subnormalities, epsilon differentiations etc.
        sigma = 1
        epsilon = 1

        return self._error(l-r,sigma=sigma,epsilon=epsilon)

    def map_product(self, expr):
        from pymbolic.primitives import Product
        assert len(expr.children) == 2

        l = self.rec(expr.children[0])
        r = self.rec(expr.children[1])

        # TODO : Define different strategies to detect subnormalities, epsilon differentiations etc.
        sigma = 1
        epsilon = 1

        return self._error(Product((l,r)),sigma=sigma,epsilon=epsilon)

    def map_quotient(self, expr):
        from pymbolic.primitives import Quotient

        l = self.rec(expr.numerator)
        r = self.rec(expr.denominator)

        # TODO : Define different strategies to detect subnormalities, epsilon differentiations etc.
        sigma = 1
        epsilon = 1

        return self._error(Quotient(l,r),sigma=sigma,epsilon=epsilon)

    def map_power(self, expr):
        from pymbolic.primitives import Power

        l = self.rec(expr.base)
        r = self.rec(expr.exponent)

        # TODO : Define different strategies to detect subnormalities, epsilon differentiations etc.
        sigma = 1
        epsilon = 1

        return self._error(Power(l,r),sigma=sigma,epsilon=epsilon)

    # TODO : Constants can also introduce rounding error
    # e.g. maybe self.rec(expr.low) and self.rex(expr.up)
    def map_interval(self, expr):
        return expr

    def map_nthroot(self, expr):
        from vodes.symbolic.expressions.nthroot import NthRoot

        l = self.rec(expr.expr)
        r = expr.n

        # TODO : Define different strategies to detect subnormalities, epsilon differentiations etc.
        sigma = 1
        epsilon = 1

        return self._error(NthRoot(l,r),sigma=sigma,epsilon=epsilon)

    def map_sin(self, expr):
        raise NotImplementedError()

    def map_cos(self, expr):
        raise NotImplementedError()

class TaylorMapper(ErrorMapper):
    def __init__(self,n:int,context:dict,min_precision:int,max_precision:int,min_exponent:int,max_exponent:int):
        from pymbolic.primitives import Variable
        from vodes.symbolic.mapper.scalar_evaluator import ScalarEvaluator
        
        super().__init__(context=context,min_precision=min_precision,max_precision=max_precision,min_exponent=min_exponent,max_exponent=max_exponent)
        
        self.err = MachineError(min_precision=min_precision, max_precision=max_precision)
        self.context[self.err.name] = Interval(self.err.bound.start, self.err.bound.end)

        self._logger = logging.getLogger(__name__)

        ## Approximations
        assert n > 0
        self.__n = n
        self.__evaluate = lambda expr: ScalarEvaluator(context=self.context)(expr)


    def round(self, te:tuple):
        from vodes.symbolic.expressions.absolute import Abs
        
        assert len(te) == 2
        (f,ts) = te

        self._logger.debug(f'Rounding : {f};{[[str(t) for t in ts]for ts in ts]}')
        # ts = [ T_1, ..., T_n, R_n \in \mathbb{M}]
        assert len(ts) == self.__n + 1

        # ROUND : (f + ts) * (1+e)
        
        # (1) TODO :Remainder via IntervalArithmetic, disable rigorousity 
        ts[-1] = [0]

        # (2) Append new expansion terms
        for i in range(self.__n - 1,-1,-1):
            # We append, as (TS_i * e^i * e)^(i+1)(0) = TS_i
            exprs = [f] if i == 0 else ts[i-1]

            ts[i].extend(
                exprs
            )

        self._logger.debug(f'-> {f};{[[str(t) for t in ts]for ts in ts]}')

        return (f,ts)

    def map_variable(self, expr):
        # TA(x) = x + 0 + .... + 0
        ts = [ [] for _ in range(0,self.__n + 1) ]

        return self.round(
            (expr,ts)
        )

    def map_constant(self, expr):
        # TA(c) = c + 0 + .... + 0
        ts = [ [] for _ in range(0,self.__n + 1) ]

        return self.round(
            (expr, ts)
        )

    def map_sum(self, expr):
        from vodes.symbolic.mapper.extended_evaluation_mapper import evaluate
        assert len(expr.children) == 2

        (fl,lts) = self.rec(expr.children[0])
        (fr,rts) = self.rec(expr.children[1])

        # TA( (F_0 + TS_0) + (F_1 + TS_1) ) = (F_0 + F_1) + (ts_00 + ts_10) + ... + (ts_0n + ts_1n)
        
        for i in range(0, self.__n + 1):
            lts[i].extend(rts[i])

        self._logger.debug(f'{expr} -> {fl + fr};{[[str(t) for t in ts]for ts in lts]}')

        res = self.round(
            (fl + fr, lts)
        )

        return res

    def map_sub(self, expr):
        from vodes.symbolic.mapper.extended_evaluation_mapper import evaluate
        assert len(expr.children) == 2

        (fl,lts) = self.rec(expr.children[0])
        (fr,rts) = self.rec(expr.children[1])

        # TA( (F_0 + TS_0) + (F_1 + TS_1) ) = (F_0 + F_1) + (ts_00 + ts_10) + ... + (ts_0n + ts_1n)
        
        for i in range(0, self.__n + 1):
            lts[i].extend([
                -rt for rt in rts[i]
            ])

        self._logger.debug(f'{expr} -> {fl - fr};{[[str(t) for t in ts]for ts in lts]}')

        return self.round(
            (fl - fr, lts)
        )

    def map_product(self, expr):
        from vodes.symbolic.mapper.extended_evaluation_mapper import evaluate
        assert len(expr.children) == 2

        (fl,lts) = self.rec(expr.children[0])
        (fr,rts) = self.rec(expr.children[1])

        # TA( (F_0 + TS_0) x (F_1 + TS_1) ) = (F_0 x F_1) + (ts_00F_0 + ts_10F_1) + ...

        # t1 = { { f_0 * t[i][j], t[i][j] \in TS_i }, TS_i \in TS_1 }
        t1 = [ [ fl * t for t in ts ] for ts in rts ]

        # t2 = { { f_1 * t[i][j], t[i][j] \in TS_i }, TS_i \in TS_0 }
        t2 = [ [ fr * t for t in ts ] for ts in lts ]

        # t3 = { { ts[i][a] * ts[j][b], i + j = z - 1, ts[i][a] \in TS_0, ts[j][b] \in TS_1}, z \in {0,...,n} }
        t3 = [ [  ] for _ in range(self.__n + 1)]
        
        for z in range(self.__n + 1):
            # 0 .. z-1
            for i in range(z):
                # z -1 .. 0
                j = (z - 1) - i

                for lt in lts[i]:
                    for rt in rts[j]:
                        t3[z].append(lt * rt)

        # t4 = Remainder
        t4 = 0 #TODO

        ts = t1
        for i in range(0, self.__n + 1):
            ts[i].extend(t2[i])
            ts[i].extend(t3[i])
        
        ts[-1].append(t4)

        return self.round(
            (fl * fr, ts)
        )
        
    def map_quotient(self, expr):
        from pymbolic.primitives import Quotient, Variable, Product

        # INV
        if isinstance(expr.numerator,int) and expr.numerator == 1:
            return self._unary_operations(
                outer=Quotient(1,Variable('x')),
                inner=self.rec(expr.denominator)
            )
        else:
            return self.rec(
                Product((expr.numerator,Quotient(1,expr.denominator)))
            )

    def map_interval(self, expr):
        raise NotImplementedError()
    ##
    # UNARY FUNCTIONS
    ##
    def _unary_operations(self, outer, inner:tuple):
        from vodes.symbolic.mapper.extended_differentiation_mapper import differentiate
        from vodes.symbolic.mapper.extended_substitution_mapper import substitute
        from vodes.symbolic.expressions.absolute import Abs

        from vodes.symbolic.mapper.interop import ExactPymbolicToSympyMapper

        (f,ts) = inner

        # First term is outer[f]
        uf = substitute(outer,{'x' : f})
        uts = [ [] for _ in range(0,self.__n + 1) ]

        derivative = outer
        for i in range(self.__n):
            derivative = differentiate(derivative,'x')
            uts[i] = [t * substitute(derivative,{'x' : f}) for t in ts[i]]
            

        # TODO : Remainder
        remainder = 0

        uts[self.__n] = [remainder]

        self._logger.debug(f'{outer};{f};{[[str(t) for t in ts]for ts in ts]} -> {uf};{[[str(t) for t in ts]for ts in uts]}')
        
        return self.round(
            (uf, uts)
        )
    def map_power(self, expr):
        from pymbolic.primitives import Power, Variable

        if not isinstance(expr.exponent, int):
            raise ValueError(f"Exponent {expr.exponent} is not supported.")
            
        return self._unary_operations(
            outer=Power(Variable('x'),expr.exponent),
            inner=self.rec(expr.base)
        )

    def map_nthroot(self, expr):
        from pymbolic.primitives import Variable
        from vodes.symbolic.expressions.nthroot import NthRoot

        if not isinstance(expr.n, int):
            raise ValueError(f"NthRoot {expr.n} is not supported.")
            
        return self._unary_operations(
            outer=NthRoot(Variable('x'),expr.n),
            inner=self.rec(expr.expr)
        )

    def map_sin(self, expr):
        from pymbolic.primitives import Variable
        from vodes.symbolic.expressions.trigonometric import sin
            
        return self._unary_operations(
            outer=sin(Variable('x')),
            inner=self.rec(expr.expr)
        )

    def map_cos(self, expr):
        from pymbolic.primitives import Variable
        from vodes.symbolic.expressions.trigonometric import cos
            
        return self._unary_operations(
            outer=cos(Variable('x')),
            inner=self.rec(expr.expr)
        )

def expand_taylor_terms(te:tuple,min_prec:int,max_prec:int,abs:bool=False):
    from functools import reduce
    assert len(te) == 2

    err = MachineError(min_precision=min_prec,max_precision=max_prec)

    inner = lambda x,i: x * Interval(err.bound.start, err.bound.end) ** (i+1)
    outer = lambda i : 1

    if abs:
        inner = lambda x,i: x
        outer = lambda i: err**(i+1)
    
    (f,ts) = te

    # Simplifies terms
    
    _error = 0

    for i in range(len(ts)):
        _error += sum(
            [
                inner(t,i) for t in ts[i]
            ]
        ) * outer(i)
    
    return (f, _error, f + _error)