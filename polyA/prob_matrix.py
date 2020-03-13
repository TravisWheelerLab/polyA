import sys
from typing import Any, Callable, Dict, List, Tuple

from .constants import NAN_STRING
from .origin_matrix import OriginMatrix

"""
A typedef that represents a probability matrix. The keys are of the form
`(i_r, i_c)` where `i_r` is the row index and `i_c` is the column index.
The value, then, is the probability of visiting that row and column.
"""
ProbMatrix = Dict[Tuple[int, int], float]


def deserialize_prob_matrix(
    matrix_lines: List[str],
) -> Tuple[ProbMatrix, OriginMatrix]:
    """
    Deserialize a probability matrix as serialized by the Perl prototype
    and Python implementations. The format of the first column ("key") is
    "<row>.<column>" (as with the support matrix). The second column
    ("prob_value") is simply a floating point probability value. The third
    column is the row in the previous column from which the value was
    derived, called "origin".
    
        key	prob_value origin
        0.0	0.0104846536699391 0
        0.1	0.0104846536699391 0
        0.2	0.0104846536699391 1
        ...

    This function assumes that the matrix was serialized with column headers.

    >>> deserialize_prob_matrix(["key\tprob_value\torigin", "1.2\t0.1\t0"])
    ({(1, 2): 0.1}, {(1, 2): 0})
    >>> deserialize_prob_matrix(["1.2\t0.1\t0"])
    Traceback (most recent call last):
      ...
    ValueError: ['1.2', '0.1', '0']
    >>> deserialize_prob_matrix(["key prob_value origin", "1.2 0.1 0"])
    ({(1, 2): 0.1}, {(1, 2): 0})
    """
    headers = matrix_lines.pop(0).split()
    if headers != ["key", "prob_value", "origin"]:
        raise ValueError(headers)

    probMatrix: ProbMatrix = {}
    originMatrix: OriginMatrix = {}

    for line in matrix_lines:
        # Ignore blank lines
        if len(line) == 0:
            continue

        [key, prob, origin] = line.split()
        [row, col] = key.split(".")
        cellPos = (int(row), int(col))
        probMatrix[cellPos] = float(prob)

        if origin == NAN_STRING:
            parsedOrigin = -1
        else:
            parsedOrigin = int(origin)
        originMatrix[cellPos] = parsedOrigin

    return (probMatrix, originMatrix)


def serialize_prob_matrix(
    prob_matrix: ProbMatrix,
    origin_matrix: OriginMatrix,
    output: Callable[[str], Any] = sys.stdout.write,
) -> None:
    """
    Serialize the probability and origin matrices, as returned by
    `fill_prob_matrix`, in a manner that is compatible with the serialization
    implemented in the Perl prototype.

    Use `output` to capture output for testing, otherwise the default value
    should be sufficient.

    >>> import io
    >>> p = {(0, 0): 1.0, (0, 1): 2.0}
    >>> o = {(0, 0): 1, (0, 1): 2}
    >>> r = io.StringIO()
    >>> serialize_prob_matrix(p, o, r.write)
    >>> print(r.getvalue())
    key prob_value origin
    0.0 1.0 1
    0.1 2.0 2
    <BLANKLINE>
    """
    output("key prob_value origin")
    for cellPos in prob_matrix.keys():
        probValue = prob_matrix[cellPos]
        originValue = origin_matrix[cellPos]

        if originValue == -1:
            originString = NAN_STRING
        else:
            originString = str(originValue)

        line = "\n%s.%s %s %s" % (
            cellPos[0],
            cellPos[1],
            probValue,
            originString,
        )
        output(line)
    output("\n")
