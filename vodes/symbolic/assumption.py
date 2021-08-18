
from typing import List

from vodes.symbolic.translations.nop_translation import NOPTranslation
from vodes.symbolic.translations.translation import Translation
from vodes.symbolic.properties.property import Property

# Custom Expression Library
from vodes.symbolic.expressions.bounded import BoundedExpression, Domain

class Assumption:
    """Class that ensures that a certain property holds and may use a given translation methods to obtain the property, if necessary.

    Args:
        property (Property): The property, that is to be given.
        translation (Translation): The translation to use, if necessary.
    """
    def __init__(self, property:Property, translation:Translation=NOPTranslation()):
        assert(property)
        assert(translation)

        self.__property = property
        self.__translation = translation

    def validate(self, exprs:List[BoundedExpression]) -> List[BoundedExpression]:
        valid = map(self.__property.verify,exprs)

        if all(valid):
            return exprs
        else:
            res = []
            for i in range(len(valid)):
                if valid[i]:
                    res.append(exprs[i])
                    continue

                tres = self.__translation(exprs[i])
                
                if all(map(self.__property.verify,tres)):
                    res.extend(tres)
                else:
                    raise ValueError(f"Neither {exprs[i]} or the result of the translation {list(map(str,tres))} fulfill the desired property.")

            return res
                