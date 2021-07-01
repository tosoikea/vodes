from sys import intern
from pymbolic.primitives import Expression, Product
from pymbolic.mapper.stringifier import PREC_NONE, StringifyMapper

class Affine(Expression):
    init_arg_names = ("noises", "values",)

    def __init__(self, noises:tuple, values:tuple):
        if len(values) - 1 != len(noises):
            raise ValueError("Invalid parameters for affine expression")

        self._noises = noises
        self._values = values

    @property
    def noises(self):
        return self._noises

    @property
    def values(self):
        return self._values

    @property
    def polynom(self):
        base = self.values[0]

        for i in range(0,len(self.noises)):
            base += self.values[i + 1] * self.noises[i]

        return base

    def __getinitargs__(self):
        return self.noises, self.values

    def make_stringifier(self, originating_stringifier=None):
        return AffineStringifyMapper()

    mapper_method = intern("map_affine")

class AffineStringifyMapper(StringifyMapper):
    def map_affine(self, expr, enclosing_prec, *args, **kwargs):
        return str(expr.polynom())