from typing import List
from vodes.symbolic import symbols

from vodes.symbolic.interval import Interval
from vodes.symbolic.symbols import Boundary, BoundedExpression, BoundedValue
from vodes.symbolic.mapper.interop import ExactPymbolicToSympyMapper
from vodes.symbolic.mapper.exact_extrema_evaluator import ExactExtremaEvaluator

# Expression library
from pymbolic.primitives import Variable, Expression

class ExactIntersectionEvaluator(ExactExtremaEvaluator):
    """Class for exactly evaluating symbolic intervals on the basis of intersections.
    TODO : Describe in depth, firstly in thesis.
    """
    def __init__(self, context: dict, symbol: Variable):
        super().__init__(context=context, symbol=symbol)

    def __intersection(self, f, g, b: Boundary):
        """Determine the real intersections of two functions for the given boundary. Intersections at the boundary values are ignored.

        Args:
            f: Function to intersect with g.
            g: Function to intersect with f.
            b: Boundary to limit the possible values for the intersections.

        Returns:
            __intersection: A (possibly empty) list of intersections between f and g.
        """
        res = self._solve(
            f=f-g,
            b=Boundary(
                lower=BoundedValue(value=b.lower.value, open=True),
                upper=BoundedValue(value=b.upper.value, open=True)
            )
        )
    
        return res

    def __intersections(self, fs:list, b: Boundary):
        """Determine the real intersections for a list of functions within the given boundary. Intersections at the boundary values are ignored.

        Args:
            fs: Functions to evaluate.
            b: Boundary to limit the possible values for the intersections.
        
        Returns:
            __intersections: A 3D array of intersections. [i][j]=[x1,x2] -> The function i (position within fs) intersects the function j at the values x1,x2 within the boundary.
        """
        intersections = [[[] for j in range(len(fs))] for i in range(len(fs))]

        for i in range(0, len(fs)):
            for j in range(0, len(fs)):
                if i == j:
                    continue

                intersections[i][j] = self.__intersection(fs[i], fs[j], b)
                intersections[j][i] = intersections[i][j]
        
        return intersections

    def __min_eval(self, x0, x1,fs:list):
        """Evaluate the minimum function within the given boundary."""
        return self.__eval(fs=fs,x0=x0, x1=x1, order=lambda a, b: a < b)

    def __max_eval(self, x0, x1,fs:list):
        """Evaluate the maximum function within the given boundary."""
        return self.__eval(fs=fs,x0=x0, x1=x1, order=lambda a, b: a > b)

    def __eval(self, order, x0, x1, fs:list):
        """Evaluate functions between two values and determine the extrema using the given order function.

        Args:
            fs: Functions to evaluate.
            x0: Start of the section
            x1: End of the section
            order: Defines the comparison between function value (e.g. <, >).
        
        Returns:
            __eval: 
        """

        values = [None] * len(fs)
        x = x0 + (x1-x0)/2

        self._logger.debug(f'Evaluating {fs} at {x} ({x0}, {x1})')

        for i in range(0, len(fs)):
            values[i] = self._boundary_eval(fs[i], BoundedValue(
                value=x,
                open=False
            ))
        
        result = (0, values[0])
        for i in range(1, len(fs)):
            if order(values[i], result[1]):
                result = (i, values[i])
            elif result[1] == values[i]:
                # As the evaluation takes places between two intersections, the functions must be equivalent, as only these cases of identical values are not contained within the intersections.
                self._logger.info(f'Encountered equivalent functions {fs[result[0]]} and {fs[i]} within the section.')

        return result[0]


    def __analysis(self, eval, b: Boundary, intersections: list, fs: list):
        """Evaluate the extremas (e.g. minima, maxima) of the given functions within the boundary and return a list of them with their respective boundary."""
        res = []

        lower = b.lower
        candidates = [i for i in range(len(fs))]
        
        # There is still domain left to analyze
        is_analyzed = False
        while not is_analyzed:
            x0 = lower.value

            # Determine the next possible intersection
            xs = [item for sublist in intersections for subsublist in sublist for item in subsublist if item > x0]
            xs.append(b.upper.value)

            x1 = min(xs)

            # Determine extremum
            extr_i = candidates[
                eval(
                    x0=x0,
                    x1=x1,
                    fs=[fs[i] for i in candidates]
                    )
                ]

            # Determine sub-domain
            upper = b.upper
            for i in range(len(intersections[extr_i])):
                xs = [item for item in intersections[extr_i][i] if item > x0 and item <= upper.value]

                if len(xs) == 0:
                    continue
                else:
                    t = min(xs)

                    if t < upper.value:
                        # may be tangent
                        candidates = [extr_i, i]

                        upper = BoundedValue(value=t,open=True)
                    else:
                        candidates.append(i)
                    
                    
            self._logger.debug(f'Extrema : {fs[extr_i]} from {lower} to {upper}')

            res.append(
                (extr_i, Boundary(lower=lower,upper=upper))
            )

            lower = BoundedValue(value=upper.value, open=not upper.open)

            # Finalize, if domain was fully analyzed
            is_analyzed = lower.value == b.upper.value
        return res

    ####
    # MIN/MAX INTERFACE
    ####
    def _get_symbolic_interval(self, pfs:list, b:Boundary):
        """Determine the upper and lower bound of the symbolic interval and return it."""
        self._logger.debug(f"Symbolic Interval Evaluation : {list(map(str,pfs))}")

        fs = [ExactPymbolicToSympyMapper()(f) for f in pfs]

        # (1) Determine intersection for possible min/max switch
        intersections = self.__intersections(fs,b)
        self._logger.debug(f'Intersections : {intersections}')

        min_res = self.__analysis(self.__min_eval, b, intersections, fs)
        self._logger.debug(f'Minima : {min_res}')
        max_res = self.__analysis(self.__max_eval, b, intersections, fs)
        self._logger.debug(f'Maxima : {max_res}')

        # Convert Bounded Values to Boundaries

        # TODO : Very inefficient, but hey
        result = []
        for minr in min_res:
            for maxr in max_res:
                combined = minr[1].intersect(maxr[1])

                if combined is None:
                    continue

                result.append(
                    BoundedExpression(
                        expression=Interval(
                            pfs[minr[0]],
                            pfs[maxr[0]]
                            ),
                        # Simplifies the expressions, however quotients are converted to floating numbers!
                        #Interval(
                        #   SympyToPymbolicMapper()(simplify(fs[minr[0]])),
                        #   SympyToPymbolicMapper()(simplify(fs[maxr[0]]))
                        #),
                        boundary=combined
                    )
                )

        self._logger.info(f"_sint: {list(map(str,pfs))}")
        self._logger.info(f"min: {min_res} , max: {max_res}")

        return result

    ####
    # SYMBOLIC INTERVAL INTERFACE
    ####
    def _sabs(self, i:Interval, b:Boundary):
        """Absolute function in the case of symbolic interval boundaries. This class uses intersections to determine the exact min/max expressions.
        
        Args:
            i (Interval): The interval to apply the absolute function on.
            b (Boundary): The boundary ('Gültigkeitsbereich') of the absolute function application.

        Returns:
            _sabs: A list of BoundedExpressions containing the result of applying the absolute function on the symbolic interval, possibly limited to subsets of the initial boundary.
        """
        lfs = self.__convert_to_positive(i.low, b)
        ufs = self.__convert_to_positive(i.up, b)

         # TODO : Improve efficiency
        result = []
        for lf in lfs:
            for uf in ufs:
                combined = lf.boundary.intersect(uf.boundary)

                if combined is None:
                    continue

                    
                result.extend(self._get_symbolic_interval(pfs=[lf.expr,uf.expr],b=combined))

        return result
        
    def _smul(self, l, r, b):
        """Interval multiplication in the case of symbolic interval boundaries. This class uses intersections to determine the exact min/max expressions.

        Args:
            l (Interval): The left parameter of the multiplication.
            r (Interval): The right parameter of the multiplication.
            b (Boundary): The boundary ('Gültigkeitsbereich') of the multiplication.

        Returns:
            _smul: A list of BoundedExpressions containing the result of the symbolic multiplication (interval), possibly limited to subsets of the initial boundary.
        """
        pfs = [
            l.low * r.low,
            l.low * r.up,
            l.up * r.low,
            l.up * r.up
        ]

        return self._get_symbolic_interval(pfs, b)

    ####
    # UTILITY FUNCTIONS
    ####
    def __convert_to_positive(self, pf:Expression, b:Boundary):
        """Determines the absolute expression.
        
        Args:
            pf (Expression): The expression to convert to an expression allowing for only positive (absolute) values.
            b (Boundary): Boundary within the expression is valid and to be converted.
            
        Returns:
            __convert_to_positive: A list of BoundedExpressions, which in sum describe the absolute expression."""
        result = []

        f = ExactPymbolicToSympyMapper()(pf)

        # Determine intersections with x axis
        xs = [b.lower]
        xs.extend(list(map(lambda x : BoundedValue(value=x, open=False),self.__intersection(f,0,b))))
        xs.append(b.upper)

        for i in range(len(xs) - 1):
            x = (xs[i].value + xs[i + 1].value) / 2
            y = self._boundary_eval(f, BoundedValue(x,False))

            bi = Boundary(
                lower = xs[i],
                upper = xs[i + 1]
            )

            if y > 0:
                result.append(
                    BoundedExpression(
                        expression=pf,
                        boundary=bi
                    )
                )
            # if y == 0 -> everywhere 0 between two points
            else:
                result.append(
                    BoundedExpression(
                        expression=(-1) * pf,
                        boundary=bi
                    )
                )

        return result