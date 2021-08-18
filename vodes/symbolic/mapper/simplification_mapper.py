import logging

# Expression Library
from pymbolic.primitives import Quotient, Sum, Product, Expression

# Expression Mapper
from pymbolic.mapper import RecursiveMapper

def simplify(expression):
    return SimplificationMapper()(expression)

class SimplificationMapper(RecursiveMapper):
    """Simple simplification mapper to collapse expressions, if possible"""
    def __init__(self):
        self._logger = logging.getLogger(__name__)

    ####
    # EXPRESSION MAPPING
    ####
    def handle_unsupported_expression(self, expr, *args, **kwargs):
        return expr

    def map_quotient(self, expr):
        return expr

    def map_constant(self, expr):
        return expr

    def map_sum(self, expr:Sum) -> Expression:
        args = list(expr.children)

        for i in range(len(args)):
            for j in range(len(args)):
                if i == j:
                    continue

                res = _simplify_add(args[i],args[j])

                if res is None:
                    continue

                args[i]=res
                args[j]=0

        args = [arg for arg in args if arg != 0]

        if len(args) == 0:
            res = 0
        elif len(args) == 1:
            res = args[0]
        else:
            res = Sum(tuple(args))

        self._logger.debug(f'{expr} -> {res}')
        return res

    def map_product(self, expr:Product) -> Expression:
        args = list(expr.children)

        for i in range(len(args)):
            for j in range(len(args)):
                if i == j:
                    continue

                res = _simplify_mul(args[i],args[j])

                if res is None:
                    continue

                args[i]=res
                args[j]=1

        args = [arg for arg in args if arg != 1]

        if len(args) == 0:
            res = 1
        elif len(args) == 1:
            res = args[0]
        else:
            res = Product(tuple(args))

        self._logger.debug(f'{expr} -> {res}')
        return res

####
# UTILITY FUNCTIONS
####
def __quotient_add(l:Quotient,r):
    res = None

    if isinstance(r,Quotient):
        res = Quotient(
            numerator=l.den * r.num + r.den * l.num,
            denominator=l.den * r.den
        )
    elif isinstance(r,int): 
        res = Quotient(
            numerator=l.num + r * l.den,
            denominator=l.den
        )

    return res

def __add_quotient_simplify(l,r):
    res = None

    if isinstance(l,Quotient):
        res = __quotient_add(l,r)
    elif isinstance(r,Quotient):
        res = __quotient_add(r,l)

    return res

def __add_natural_simplify(l,r):
    res = None

    if isinstance(l,int) and isinstance(r,int):
        res = l + r
    
    return res

def _simplify_add(l,r):
    _simplifications = [
        __add_natural_simplify,
        __add_quotient_simplify
    ]

    res = None

    for _simplification in _simplifications:
        res = _simplification(l,r)

        if not(res is None):
            break

    return res

def __quotient_mul(l:Quotient,r):
    res = None
    
    if isinstance(r,Quotient):
        res = Quotient(
            numerator=l.num * r.num,
            denominator=l.den * r.den
        )
    elif isinstance(r,int): 
        res = Quotient(
            numerator=l.num * r,
            denominator=l.den
        )

    return res

def __mul_quotient_simplify(l,r):
    res = None

    if isinstance(l,Quotient):
        res = __quotient_mul(l,r)
    elif isinstance(r,Quotient):
        res = __quotient_mul(r,l)

    return res

def __mul_natural_simplify(l,r):
    res = None

    if isinstance(l,int) and isinstance(r,int):
        res = l * r
    
    return res

def _simplify_mul(l,r):
    _simplifications = [
        __mul_natural_simplify,
        __mul_quotient_simplify
    ]

    res = None

    for _simplification in _simplifications:
        res = _simplification(l,r)

        if not(res is None):
            break

    return res