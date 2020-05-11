import sys
import numpy as np

from math import inf
from typing import Any, Callable, Dict, List, Tuple

from .conf_score_matrix import ConfScoreMatrix

SupportMatrix = Dict[Tuple[int, int], float]
"""
Hash implementation of sparse 2D matrix used in pre-DP calculations. Key is
tuple[int, int] that maps row, col with the value held in that cell of matrix. Rows are 
subfamilies in the input alignment file, cols are nucleotide positions in the alignment.
Each cell in matrix is the support score (or average confidence value) for the following 
chunksize cells in the confidence matrix.
"""


CollapsedSupportMatrix = Dict[str, float]
"""
One the support matrix is collapsed it becomes a mapping from
subfamily name into the total support value for that subfamily.

TODO (Kaitlin): Confirm this and feel free to edit
"""


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


def fill_support_matrix(
    conf_matrix: ConfScoreMatrix, n_rows: int, non_empty_columns: List[int]
) -> SupportMatrix:
    support_matrix: SupportMatrix = {}

    """
    Fills support score matrix using values in conf matrix. Score for subfam row 
    at position col is sum of all confidences for subfam row for that column and 
    the following chunksize-1 columns - normalized by dividing by number of segments
    
    Ex: support_matrix[0,0] is sum of conf_matrix[0,0] to conf_matrix[0,30], divided by 31
    """
    for i in range(n_rows):
        for col in range(len(non_empty_columns)):
            j = non_empty_columns[col]

            if (i, j) in conf_matrix:

                num = 0
                summ = 0
                numsegments = 0
                while num < chunksize:
                    if (i, j + num) in conf_matrix:
                        summ += conf_matrix[i, j + num]
                        numsegments += 1

                    num += 1

                support_matrix[i, j] = summ / numsegments

    return support_matrix


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
