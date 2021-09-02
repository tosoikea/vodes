from abc import ABC, abstractmethod
from typing import List
from vodes.symbolic.mapper.interop import ExactPymbolicToSympyMapper
from vodes.symbolic.expressions.bounded import BoundedExpression, MachineError

class Solution(ABC):
    """Class to represent solution to an error analysis problem"""
    @abstractmethod
    def error(self, prec:int):
        pass

    @property
    def name(self):
        """Get or set the name of the solution"""
        pass

class PseudoExactSolution(Solution):
    def __calculate(self,prec):
        from mpmath import mp
        t = mp.prec
        try:
            mp.prec = prec
            return self.func()
        finally:
            mp.prec = t

    def __init__(self, func, name:str="Exact"):
        self.func = func
        # 113 ~ exact solution
        self.expected = self.__calculate(113)
        self.__name = name

    @property
    def name(self):
        """Get the name of the solution"""
        return self.__name

    def error(self, prec:int):
        return abs(self.expected -  self.__calculate(prec))

class PseudoExactIntervalSolution(Solution):
    def __calculate(self,prec):
        from mpmath import iv
        t = iv.prec
        try:
            iv.prec = prec
            return self.func()
        finally:
            iv.prec = t

    def __init__(self, func, name:str="Exact"):
        self.func = func
        # 113 ~ exact solution
        self.expected = self.__calculate(113)
        self.__name = name

    @property
    def name(self):
        """Get the name of the solution"""
        return self.__name

    def error(self, prec:int):
        from mpmath import mpf

        iv = abs(self.expected - self.__calculate(prec))

        return mpf(max([
            iv.a,
            iv.b
        ])._mpi_[0])

class AnalysisSolution(Solution):
    def __init__(self, bexprs:List[BoundedExpression], name:str):
        super().__init__()

        self.bexprs = [
            BoundedExpression(
                expression=ExactPymbolicToSympyMapper()(bexpr.expr),
                boundary=bexpr.bound
            )for bexpr in bexprs]

        self.__name = name

    @property
    def name(self):
        """Get the name of the solution"""
        return self.__name

    def error(self, prec: int):
        from sympy import Float
        from vodes.symbolic.mapper.extended_evaluation_mapper import evaluate
        precision = 256

        err = evaluate(MachineError(min_precision=prec,max_precision=prec).bound.end)

        for bexpr in self.bexprs:
            if not bexpr.bound.contains(err):
                continue
            
            if len(bexpr.expr.free_symbols) == 0:
                return bexpr.expr.evalf(precision)
            elif len(bexpr.expr.free_symbols) == 1:
                sym = list(bexpr.expr.free_symbols)[0]
                return bexpr.expr.subs(sym,Float(err,precision)).evalf(precision)
            else:
                raise ValueError("Encountered too many free variables.")

def show(solutions:List[Solution],min_prec:int=11,max_prec:int=53):
    from matplotlib import pyplot
    assert max_prec > min_prec
    
    xs = range(min_prec,max_prec+1)

    for sol in solutions:
        ys = []

        for x in xs:
            ys.append(sol.error(prec=x))

        pyplot.plot(xs,ys, label=sol.name)

    pyplot.grid(True)
    pyplot.legend()
    pyplot.show()