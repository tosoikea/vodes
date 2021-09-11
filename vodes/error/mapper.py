from functools import reduce
import logging

from vodes.symbolic.mapper.interop import ExactPymbolicToSympyMapper, ExactSympyToPymbolicMapper
from vodes.symbolic.expressions.nthroot import NthRoot
from vodes.symbolic.expressions.primitives import Subtraction
from vodes.symbolic.expressions.trigonometric import sin, cos
from vodes.symbolic.expressions.bounded import DummyVariable, MachineError, SmallestMachineNumber
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
    def __init__(self,context:dict,min_precision:int,max_precision:int,min_exponent:int,max_exponent:int):
        from pymbolic.primitives import Variable
        from vodes.symbolic.mapper.taylor_evaluator import TaylorEvaluator
        
        super().__init__(context=context,min_precision=min_precision,max_precision=max_precision,min_exponent=min_exponent,max_exponent=max_exponent)

        self._logger = logging.getLogger(__name__)

        ## Approximations

    ###
    # TAYLOR FORM INTERFACE
    ###
    def __compact(self, expr):
        # TODO : Provide proper compaction within this framework or implement it within pymbolic
        from vodes.symbolic.mapper.interop import ExactPymbolicToSympyMapper, ExactSympyToPymbolicMapper
        
        return ExactSympyToPymbolicMapper()(
            ExactPymbolicToSympyMapper()(
                expr
            ).doit()
        )

    def _bound(self,expr):
        import copy
        from vodes.symbolic.expressions.absolute import Abs
        from vodes.symbolic.mapper.scalar_evaluator import evaluate
        ## TODO : Improve the performance of the symbolic interval evaluator to allow usage within the bounding process
        ##        This in turn allows for the machine epsilon variables to be kept and allows for better approximation over ranges of precisions
        _context = copy.deepcopy(self.context)
        _context[self.eps.name] = Interval(self.eps.bound.start,self.eps.bound.end)

        res = evaluate(
            Abs(expr),
            context=_context
        )

        if len(res) > 1:
            raise ValueError("Unexpected divergence of expressions during bounding")

        return res[0].expr

    def _tadd(self, sf:tuple, tf:tuple):
        """Computes the addition of two taylor forms"""
        (s,(b1,b2)) = sf
        (t,(c1,c2)) = tf

        s_k = len(b1)

        op_f = self.__compact(s + t)

        op_c1 = b1
        op_c1.extend(c1)

        op_c2 = [
            [
                0 for _ in range(len(op_c1))
            ] 
            for _ in range(len(op_c1))
        ]

        self._logger.debug(f'{s};{list(map(str,b1))};{[list(map(str,cs)) for cs in b2]}')
        self._logger.debug(f'{t};{list(map(str,c1))};{[list(map(str,cs)) for cs in c2]}')

        for i in range(len(op_c2)):
            for j in range(len(op_c2)):

                if i < s_k and j < s_k:
                    op_c2[i][j] = b2[i][j]
                elif i >= s_k and j >= s_k:
                    op_c2[i][j] = c2[i - s_k][j - s_k]

        # Ensure correctness
        assert len(op_c2[-1]) == len(op_c2[0]) == len(op_c1)
        return (op_f, (op_c1,op_c2))

    def _tsub(self, sf:tuple, tf:tuple):
        """Computes the subtraction of two taylor forms"""
        
        (s,(b1,b2)) = sf
        (t,(c1,c2)) = tf

        s_k = len(b1)

        op_f = self.__compact(s - t)

        op_c1 = b1
        op_c1.extend([-c for c in c1])

        op_c2 = [
            [
                0 for _ in range(len(op_c1))
            ] 
            for _ in range(len(op_c1))
        ]

        for i in range(len(op_c2)):
            for j in range(len(op_c2)):

                if i < s_k and j < s_k:
                    op_c2[i][j] = b2[i][j]
                elif i >= s_k and j >= s_k:
                    op_c2[i][j] = -c2[i - s_k][j - s_k]

        # Ensure correctness
        assert len(op_c2[-1]) == len(op_c2[0]) == len(op_c1)
        return (op_f, (op_c1,op_c2))

    def _tmul(self, sf:tuple, tf:tuple):
        """Computes the multiplication of two taylor forms"""
        from vodes.symbolic.expressions.absolute import Abs
        from vodes.symbolic.mapper.comparison_evaluator import evaluate

        (s,(b1,b2)) = sf
        (t,(c1,c2)) = tf

        s_k = len(b1)
        t_k = len(c1)
        c_k = s_k + t_k

        mul_f = self.__compact(s * t)

        # (1). Extend first order 
        mul_c1 = [0 for _ in range(s_k + t_k)]
        for i in range(c_k):
            if i < s_k:
                mul_c1[i] = self.__compact(t * b1[i])
            else:
                mul_c1[i] = self.__compact(s * c1[i - s_k])

        # (2). Extend second order
        mul_c2 = [
            [
                0 for _ in range(c_k)
            ] 
            for _ in range(c_k)
        ]

        for i in range(c_k):
            for j in range(c_k):

                if i < s_k:
                    if j < s_k:
                        mul_c2[i][j] = self.__compact(t * b2[i][j])
                    else:
                        mul_c2[i][j] = self.__compact(b1[i] * c1[j - s_k])
                else:
                    if j < s_k:
                        mul_c2[i][j] = self.__compact(c1[i - s_k] * b1[j])
                    else:
                        mul_c2[i][j] = self.__compact(s * c2[i - s_k][j - s_k])

        # (3). Bound Remainder
        _em1 = 0
        _em2 = 0
        _em3 = 0

        for i in range(s_k):
            # M1
            for m in range(t_k):
                for g in range(t_k):
                    _em1 += Abs(b1[i] * c2[m][g])

            for j in range(s_k):
                for m in range(t_k):
                    # M2
                    _em1 += Abs(b2[i][j] * c1[m])

                    # M3
                    for g in range(t_k):
                        _em3 += Abs(b2[i][j] * c2[m][g])

        m1 = self._bound(_em1 + _em2)
        m3 = self._bound(_em3)

        # (4). Add with remainder form
        (op_f,(op_c1,op_c2)) = self._tadd(
            (mul_f,(mul_c1,mul_c2)),
            (0,([
                    self.__compact(
                        Power(self.eps,2) * m1 + Power(self.eps,3) * m3
                    )
                ],[[0]]))
        )        
        
        self._logger.debug(f'-> {op_f};{list(map(str,op_c1))};{[list(map(str,cs)) for cs in op_c2]}')
        return (op_f,(op_c1,op_c2))

    def _tround(self, te:tuple):
        """Computes the taylor form for the rounding of an existing taylor form"""
        from vodes.symbolic.expressions.absolute import Abs
        
        assert len(te) == 2
        (f,(c1,c2)) = te

        self._logger.debug(f'Rounding : {f};{list(map(str,c1))};{[list(map(str,cs)) for cs in c2]}')

        # ts : [C_0,C_1]
        assert not (f is None)

        # ROUND : (f + ts) * (1+e)

        # (1). Bound remainder
        m = self._bound(
            sum(
                [sum([Abs(c) for c in cs]) for cs in c2]
            ),
        )

        # (2). Append e_{k+1} for function
        c2.append([])
        for i in range(len(c2) - 1):
            c2[i].append(c1[i])
            c2[-1].append(c1[i])

        c2[-1].append(0)
        c1.append(f)

        # Ensure correctness
        assert len(c2[-1]) == len(c2[0]) == len(c1)

        # (3). Add with remainder form
        (op_f,(op_c1,op_c2)) = self._tadd(
            (f,(c1,c2)),
            (0,([
                self.__compact(Power(self.eps,2) * m)
                ],[[0]]))
        )

        self._logger.debug(f'-> {op_f};{list(map(str,op_c1))};{[list(map(str,cs)) for cs in op_c2]}')
        return (op_f,(op_c1,op_c2))

    def _tfunc(self, te:tuple, g, g_1, g_2, g_3):
        from pymbolic.primitives import Quotient
        from vodes.symbolic.mapper.extended_substitution_mapper import substitute
        from vodes.symbolic.expressions.absolute import Abs

        (t,(c1,c2)) = te
        self._logger.debug(f'Applying {g} to {t};{list(map(str,c1))};{[list(map(str,cs)) for cs in c2]}')

        t_k = len(c1)

        g_f = self.__compact(substitute(g,{'x':t}))

        # (1). Bound remainder
        _em0 = sum([Abs(c) for c in c1])
        _em1 = sum([
            sum([Abs(c2[i][j] + Abs(c2[j][i])) for j in range(t_k)]) for i in range(t_k)
        ])

        (_,_,h_t) = expand_taylor_terms(te,min_prec=self.min_precision,max_prec=self.max_precision,abs=False)
        _em2 = substitute(g_2,{'x':h_t})
        _em3 = substitute(g_3,{'x':h_t})
        
        m0 = self._bound(_em0)
        m1 = self._bound(_em1)
        m2 = self._bound(_em2)
        m3 = self._bound(_em3)

        m_a = 2 * t_k * m1 * m2
        m_b = t_k**2 * m0 * m3
        m_c = t_k**2 * m1 * m3

        # (2). Extend first order 
        for i in range(t_k):
            c1[i] = self.__compact(substitute(g_1,{'x':t}) * c1[i])

        # (3). Extend second order
        for i in range(t_k):
            for j in range(t_k):
                c2[i][j] = self.__compact(substitute(g_2,{'x':t}) * c1[i] + substitute(g_1,{'x':t}) * (c2[i][j] + c2[j][i]))

        # (4). Add with remainder form
        (op_f,(op_c1,op_c2)) = self._tadd(
            (g_f,(c1,c2)),
            (0,((
                [
                    # 1/6 * (e^2M_a + e^2M_b + e^3M_c)
                    self.__compact(
                        Quotient(1,6) * 
                            (Power(self.eps,2) * m_a 
                                + Power(self.eps,2) * m_b 
                                + Power(self.eps,3) * m_c)
                    )
                ]
            ),[[0]]))
        )

        self._logger.debug(f'-> {op_f};{list(map(str,op_c1))};{[list(map(str,cs)) for cs in op_c2]}')
        return (op_f,(op_c1,op_c2))
        
    def _unary_operations(self, tf:tuple, g):
        from vodes.symbolic.mapper.extended_differentiation_mapper import differentiate
        from vodes.symbolic.mapper.extended_substitution_mapper import substitute

        # (1) Determine derivatives
        g_1 = differentiate(g,'x')
        g_2 = differentiate(g_1,'x')
        g_3 = differentiate(g_2,'x')

        return self._tround(
            self._tfunc(
                tf,
                g,
                g_1,
                g_2,
                g_3
            )
        )

    ###
    # UTILITY FUNCTIONS
    ###
    def is_power_two(self, expr):
        #TODO : Negative exponents e.g. x=2^-5?

        if isinstance(expr,int):
            #is power of two? (src : http://graphics.stanford.edu/~seander/bithacks.html#DetermineIfPowerOf2)
            v = abs(expr)
            return v and not(v & (v-1))
        else:
            return False

    ###
    # EXPRESSION MAPPING
    ###
    def map_variable(self, expr):
        # TA(x) = x
        ts = ([],[])

        return self._tround(
            (expr,ts)
        )

    def map_constant(self, expr):
        # Improved error model
        b = (expr,([],[]))

        if self.is_power_two(expr):
            return b
        else:
            return self._tround(
                b
            )

    def map_sum(self, expr):
        from vodes.symbolic.mapper.extended_evaluation_mapper import evaluate
        assert len(expr.children) == 2

        return self._tround(
            self._tadd(
                self.rec(expr.children[0]),
                self.rec(expr.children[1])
            )
        )
        
    def map_sub(self, expr):
        from vodes.symbolic.mapper.extended_evaluation_mapper import evaluate
        assert len(expr.children) == 2

        return self._tround(
            self._tsub(
                self.rec(expr.children[0]),
                self.rec(expr.children[1])
            )
        )

    def map_product(self, expr):
        from vodes.symbolic.mapper.extended_evaluation_mapper import evaluate
        assert len(expr.children) == 2

        return self._tround(
            self._tmul(
                self.rec(expr.children[0]),
                self.rec(expr.children[1])
            )
        )

    def map_interval(self, expr):
        raise NotImplementedError()
    
    def map_power(self, expr):
        from pymbolic.primitives import Power, Variable

        if not isinstance(expr.exponent, int):
            raise ValueError(f"Exponent {expr.exponent} is not supported.")
            
        return self._unary_operations(
            tf  =   self.rec(expr.base),
            g   =   Power(Variable('x'),expr.exponent)
        )

    ##
    # UNARY EXPRESSIONS
    ##
    def map_quotient(self, expr):
        from pymbolic.primitives import Quotient, Variable, Product

        # INV
        if isinstance(expr.numerator,int) and expr.numerator == 1:
            return self._unary_operations(
                tf  =   self.rec(expr.denominator),
                g   =   Quotient(1,Variable('x'))
            )
        else:
            return self.rec(
                Product((expr.numerator,Quotient(1,expr.denominator)))
            )

    def map_nthroot(self, expr):
        from pymbolic.primitives import Variable
        from vodes.symbolic.expressions.nthroot import NthRoot

        if not isinstance(expr.n, int):
            raise ValueError(f"NthRoot {expr.n} is not supported.")
            
        return self._unary_operations(
            tf  =   self.rec(expr.expr),
            g   =   NthRoot(Variable('x'),expr.n)
        )

    def map_sin(self, expr):
        from pymbolic.primitives import Variable
        from vodes.symbolic.expressions.trigonometric import sin

        return self._unary_operations(
            tf  =   self.rec(expr.expr),
            g   =   sin(Variable('x'))
        )

    def map_cos(self, expr):
        from pymbolic.primitives import Variable
        from vodes.symbolic.expressions.trigonometric import cos
            
        return self._unary_operations(
            tf  =   self.rec(expr.expr),
            g   =   cos(Variable('x'))
        )


def expand_taylor_terms(te:tuple,min_prec:int,max_prec:int,abs:bool=False,separated:bool=False):
    """Converts the taylor form to its written out form"""
    from vodes.symbolic.expressions.absolute import Abs
    
    err = MachineError(min_precision=min_prec,max_precision=max_prec)

    assert len(te) == 2
    (f,(c1,c2)) = te

    inner = lambda x,i: x * Interval(-err, err) ** i
    outer = lambda x,i: x

    if abs:
        inner = lambda x,i: Abs(x)
        outer = lambda x,i: x * err**i
    
    # (1) Order one
    _error_one = sum([
            inner(c,1) for c in c1
        ])

    # (2) Order two
    _error_two = 0
    for i in range(len(c1)):
        # (i,j) (j,i) have same error variables
        for j in range(i,len(c1)):
            _error_two += inner(c2[i][j] + c2[j][i],2)
    
   
    if separated:
        return (f,
            [(_error_one,outer(1,1)),(_error_two,outer(1,2))]
        )
    else:
        _error = outer(_error_one,1) + outer(_error_two,2)
        return (f, _error, f + _error)