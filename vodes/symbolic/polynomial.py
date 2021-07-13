from abc import abstractmethod, ABC
from sys import intern

from pymbolic.mapper.stringifier import StringifyMapper
from vodes.symbolic.symbols import MachineError
from numpy.lib.arraysetops import isin
from pymbolic.polynomial import Polynomial
from pymbolic.primitives import Quotient, quotient

from pymbolic.interop.sympy import SympyToPymbolicMapper, PymbolicToSympyMapper 
from sympy import symbols, AccumBounds

class ApproxPolynomial(ABC,Polynomial):
    mapper_method = intern("map_approx_polynomial")

    def make_stringifier(self, originating_stringifier=None):
        return ApproxPolynomialStringifyMapper()

    def __truediv__(self, other):
        if not isinstance(other, Polynomial):
            return 1/other * self

        po = Polynomial(other.Base, other.Data, other.unit, other.VarLess)
        ps = Polynomial(self.Base, self.Data, self.unit, self.VarLess)

        q,r = divmod(ps,po)

        # no remainder
        if r.degree == -1:
            print("T1")
            return MachinePolynomial(q.Base,q.Data)
        # approximation of remainder using laurent series
        else:
            print("T2")
            return MachinePolynomial(q.Base,q.Data) + self._approx(r,po)

    @abstractmethod
    def _approx_remainder(self, r, d):
        pass

class MachinePolynomial(ApproxPolynomial):
    """
    Knowledge : free variable is MachineVariable
    """
    def _approx_remainder(self, r, d):
        print("T3")
        if not isinstance(r.Base, MachineError):
            raise TypeError(f'Cannot approximate base {r.Base}')

        expr = PymbolicToSympyMapper()(Quotient(r,d))
        print(expr)
        # TODO : get machineError -> sympy only once
        bounds = expr.subs(PymbolicToSympyMapper()(MachineError()), AccumBounds(0,1))
        print(bounds)
        return SympyToPymbolicMapper()(bounds.max)


class ApproxPolynomialStringifyMapper(StringifyMapper):
    def map_approx_polynomial(self, expr, enclosing_prec, *args, **kwargs):
        return str(Polynomial(expr.Base, expr.Data))