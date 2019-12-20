import sys
from typing import Any, Callable

from .origin_matrix import OriginMatrix
from .prob_matrix import ProbMatrix


def serialize_prob_matrix(
    probMatrix: ProbMatrix,
    originMatrix: OriginMatrix,
    output: Callable[[str], Any] = sys.stdout.write,
) -> None:
    """
    Serialize the probability and origin matrices, as returned by
    `fill_prob_matrix`, in a manner that is compatible with the serialization
    implemented in the Perl prototype.

    >>> import io
    >>> p = {(0, 0): 1.0, (0, 1): 2.0}
    >>> o = {(0, 0): 1, (0, 1): 2}
    >>> r = io.StringIO()
    >>> serialize_prob_matrix(p, o, r.write)
    >>> print(r.getvalue())
    key prob_value origin
    0.0 1.0 1
    0.1 2.0 2
    """
    output('key prob_value origin')
    for cellPos in probMatrix.keys():
        probValue = probMatrix[cellPos]
        originValue = originMatrix[cellPos]
        line = '\n%s.%s %s %s' % (cellPos[0], cellPos[1], probValue, originValue)
        output(line)
