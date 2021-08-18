from sys import intern

class Infinity:
    def __str__(self) -> str:
        return f'oo'

    def __eq__(self, other):
        return isinstance(other, Infinity)

    def __ne__(self, other):
        return not self == other

    def __lt__(self, _):
        return False

    def __le__(self, other):
        if isinstance(other, Infinity):
            return True
        else:
            return False

    def __gt__(self, other):
        if isinstance(other, Infinity):
            return False
        else:
            return True

    def __ge__(self, _):
        return True

    def __hash__(self):
        return hash(Infinity)

class NegativeInfinity:
    def __str__(self) -> str:
        return f'-oo'

    def __eq__(self, other):
        return isinstance(other, NegativeInfinity)

    def __ne__(self, other):
        return not self == other

    def __lt__(self, other):
        if isinstance(other, NegativeInfinity):
            return False
        else:
            return True

    def __le__(self, _):
        return True

    def __gt__(self, _):
        return False

    def __ge__(self, other):
        if isinstance(other, NegativeInfinity):
            return True
        else:
            return False

    def __hash__(self):
        return hash(NegativeInfinity)