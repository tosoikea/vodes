from sympy import symbols, shape
from sympy.tensor.functions import NoShapeError

class Problem:
    ##
    # The problem has to be a symbolic expression with the following restrictions.
    # y'(t) = 3 * y(t) -> Problem(3 * symbols('y'))
    # u'(t) = 3 * u(t) -> Problem(3 * symbols('y'))
    # y'(t) = -3 * y(t) + t ** y(t) -> Problem(-3 * symbols('y') + symbols('t') ** symbols('y'))
    ##
    def __init__(self, problem, iv, interval):
        self.__problem = problem
        self.__iv = iv
        self.__interval = interval

        # Differentiate between system of equations and single equation (of first order)
        try:
            (_,_) = shape(self.get_expression(*symbols('y t')))
            self.matrix = True
        except NoShapeError:
            self.matrix = False
    
    def get_expression(self, y, t):
        return self.__problem.subs({
            symbols('y') : y,
            symbols('t') : t
        })

    def start(self):
        return self.__interval[0]

    def end(self):
        return self.__interval[1]

    def initial (self):
        return self.__iv

    def validate(self):
        if len(self.__interval) != 2:
            exit(-1)