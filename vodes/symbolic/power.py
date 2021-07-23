from pymbolic.primitives import Power as PE

class Power(PE):
    def __compare_power(self, other, f):
        self_val = self.base ** self.exponent

        if isinstance(other,PE):
            # a) same base => ordering via exponent
            if self.base == other.base:
                return f(self.exponent, other.exponent)
            else:
            # b) differentiating base
                return f(self_val, other.base ** other.exponent)
        else:
            return f(self_val, other)

    def __eq__(self, other):
        if isinstance(other, PE):
            return self.base == other.base and self.exponent == other.exponent

        return False

    def __ne__(self, other):
        return not self == other

    def __lt__(self, other):
        return self.__compare_power(other, f=lambda a,b : a < b)

    def __le__(self, other):
        return self.__compare_power(other, f=lambda a,b : a <= b)

    def __gt__(self, other):
        return self.__compare_power(other, f=lambda a,b : a > b)

    def __ge__(self, other):
        return self.__compare_power(other, f=lambda a,b : a >= b)

    def __hash__(self):
        return hash((self.base, self.exponent))