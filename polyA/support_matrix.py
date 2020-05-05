import sys
import numpy as np

from math import inf
from typing import Any, Callable, Dict, List, Tuple

"""
Typedef to represent a sparse support matrix implemented as a
dictionary (for now).
"""
SupportMatrix = Dict[Tuple[int, int], float]


def deserialize_support_matrix(matrix_lines: List[str]) -> SupportMatrix:
    """
    Deserialize a support matrix as serialized by the Perl prototype
    implementation. The format of the first column ("key") is "<row>.<column>".
    The second column ("support_value") is simply a floating point value.
    
        key	support_value
        0.0	0.0104846536699391
        0.1	0.0104846536699391
        0.2	0.0104846536699391
        ...

    This function assumes that the matrix was serialized with column headers.

    >>> deserialize_support_matrix(["key\tsupport_value", "1.2\t0.1"])
    {(1, 2): 0.1}
    >>> deserialize_support_matrix(["1.2\t0.1"])
    Traceback (most recent call last):
      ...
    ValueError: ['1.2', '0.1']
    >>> deserialize_support_matrix(["key support_value", "1.2 0.1"])
    {(1, 2): 0.1}
    """
    # Verify the headers
    headers = matrix_lines.pop(0).split()
    if headers != ["key", "support_value"]:
        raise ValueError(headers)

    support_matrix: SupportMatrix = {}

    # Populate the actual matrix values
    for line in matrix_lines:
        # Skip empty lines
        if len(line) == 0:
            continue

        [key, value] = line.split()
        [row, col] = key.split(".")
        cell_pos = (int(row), int(col))
        support_matrix[cell_pos] = float(value)

    return support_matrix


def fill_support_matrix() -> SupportMatrix:
    pass


def serialize_support_matrix(
    support_matrix: SupportMatrix,
    output: Callable[[str], Any] = sys.stdout.write,
) -> None:
    """
    Serialize the support matrix in a manner that is compatible with the
    serialization implemented in the Perl prototype.

    Use `output` to capture output for testing, otherwise the default value
    should be sufficient.

    >>> import io
    >>> s = {(0, 0): 1.0, (0, 1): 2.0}
    >>> r = io.StringIO()
    >>> serialize_support_matrix(s, r.write)
    >>> print(r.getvalue())
    key support_value
    0.0 1.0
    0.1 2.0
    <BLANKLINE>
    """
    output("key support_value")
    for cellPos in support_matrix:
        support_value = support_matrix[cellPos]
        line = "\n%s.%s %s" % (cellPos[0], cellPos[1], support_value)
        output(line)
    output("\n")


def support_matrix_dims(
    support_matrix: SupportMatrix,
) -> Tuple[int, int, float, float]:
    """
    Returns a 4-tuple of the number of rows and columns in the matrix,
    respectively, followed by the minimum and maximum values in the
    matrix, again respectively.
    """
    rows: int = 0
    cols: int = 0
    min_value: float = inf
    max_value: float = -inf
    for (row, col), value in support_matrix.items():
        if row >= rows:
            rows = row
        if col >= cols:
            cols = col
        if value < min_value:
            min_value = value
        if value > max_value:
            max_value = value
    return rows, cols, min_value, max_value
