from abc import ABC, abstractmethod
from typing import List
from vodes.symbolic.mapper.interop import ExactPymbolicToSympyMapper
from vodes.symbolic.expressions.bounded import BoundedExpression, MachineError

class Solution(ABC):
    """Class to represent solution to an error analysis problem"""

    @abstractmethod
    def error(self, prec:int):
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

    def __init__(self, func):
        self.func = func
        # 113 ~ exact solution
        self.expected = self.__calculate(113)

    def error(self, prec:int):
        return abs(self.expected -  self.__calculate(prec))

class AnalysisSolution(Solution):
    def __init__(self, bexprs:List[BoundedExpression]):
        super().__init__()

        self.bexprs = [
            BoundedExpression(
                expression=ExactPymbolicToSympyMapper()(bexpr.expr),
                boundary=bexpr.bound
            )for bexpr in bexprs]

    def error(self, prec: int):
        from sympy import Float
        from vodes.symbolic.mapper.extended_evaluation_mapper import evaluate

        err = evaluate(MachineError(min_precision=prec,max_precision=prec).bound.end)

        for bexpr in self.bexprs:
            if not bexpr.bound.contains(err):
                continue

            if len(bexpr.expr.free_symbols) == 0:
                return bexpr.expr.evalf(512)
            elif len(bexpr.expr.free_symbols) == 1:
                sym = list(bexpr.expr.free_symbols)[0]
                return bexpr.expr.subs(sym,Float(err,512)).evalf(512)
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
    
        pyplot.plot(xs,ys)

    pyplot.grid(True)
    pyplot.show()