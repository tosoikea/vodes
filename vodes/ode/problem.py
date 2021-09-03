from vodes.symbolic.expressions.matrix import Matrix

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
        self.matrix = isinstance(iv,Matrix)

        if self.matrix:
            raise ValueError("Matrix is not yet supported.")
    
    def get_expression(self, y, t):
        from vodes.symbolic.mapper.extended_substitution_mapper import substitute

        return substitute(self.__problem, variable_assignments={
            'y' : y,
            't' : t
        })

    @property
    def start(self):
        return self.__interval[0]

    @property
    def end(self):
        return self.__interval[1]

    @property
    def initial (self):
        return self.__iv

    def validate(self):
        assert len(self.__interval) == 2
        assert self.__interval[0] < self.__interval[1]