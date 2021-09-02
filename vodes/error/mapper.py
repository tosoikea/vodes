from vodes.symbolic.expressions.primitives import Subtraction
from vodes.symbolic.expressions.trigonometric import sin, cos
from vodes.symbolic.expressions.nthroot import NthRoot
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

class AffineMapper(ErrorMapper):
    def __init__(self,context:dict,min_precision:int,max_precision:int,min_exponent:int,max_exponent:int):
        super().__init__(context=context,min_precision=min_precision,max_precision=max_precision,min_exponent=min_exponent,max_exponent=max_exponent)
    
        self.__id = 0
        self.__noises = []
        self.__exact = {
            "sigma" : 0
        }

        self.__ids = {}

    @property
    def noises(self):
        """Get the noise variables that were created using this mapper"""
        return self.__noises

    @property
    def exact(self):
        """Get the substitutions of all created noise variables to zero"""
        return self.__exact

    @property
    def ids(self):
        """Get the mapping from ids to expressions"""
        return self.__ids

    def _error(self, expr, sigma:int, epsilon:int):
        expr_id = str(expr)

        if not expr_id in self.__ids:
            # TODO : This would need to be tweaked in case of parallelization
            t = self.__id
            self.__id = self.__id + 1
            self.__ids[expr_id] = t

        eps = Variable(name=f'eps_{self.__ids[expr_id]}')

        if not eps.name in self.__exact:
            self.__noises.append(eps)
            self.__exact[eps.name] = 0

        res = expr * (1 + epsilon * eps)

        # TODO : Underflow
        # + sigma * self.sigm.bound.end
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
    def map_interval(self, expr):
        return expr

    def map_nthroot(self, expr):
        raise NotImplementedError()

    def map_sin(self, expr):
        raise NotImplementedError()

    def map_cos(self, expr):
        raise NotImplementedError()