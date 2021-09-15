from typing import Dict, List, Tuple

from .confidence_cm import confidence_cm
from .matrices import ConfidenceMatrix
from .performance import timeit


@timeit()
def fill_confidence_matrix(
    column_count: int,
    subfam_counts: Dict[str, float],
    subfams: List[str],
    active_cells: Dict[int, List[int]],
    align_matrix: Dict[Tuple[int, int], float],
) -> ConfidenceMatrix:
    """
    Fills confidence matrix from alignment matrix. Each column in the alignment
    matrix is a group of competing annotations that are input into
    confidence_cm, the output confidence values are used to populate
    confidence_matrix.

    Inputs:

    columns - non empty matrix column indices
    subfam_counts - mapping of subfamily names to their prior counts
    subfams - list of subfamily names taken from the original alignments
    active_cells - map of column indices (from columns) into list of non-empty
    rows for the given column
    align_matrix - the alignment matrix from fill_align_matrix

    Outputs:

    confidence_matrix - hash implementation of sparse 2D matrix used in pre-DP
    calculations. Key is tuple[int, int] that maps row, col with the value held
    in that cell of matrix. Rows are subfamilies in the input alignment file,
    cols are nucleotide positions in the alignment. Each cell in matrix is the
    confidence score calculated from all the alignment scores in a column of the
    alignment matrix.

    >>> align_mat = {
    ...     (0, 0): 10, (0, 1): 10, (0, 2): 10, (0, 3): 10, (0, 4): 10,
    ...     (1, 1): 42, (1, 2): 41, (1, 3): 41, (2, 1): 45, (2, 2): 43, (2, 3): 39,
    ... }
    >>> active = {0: [0], 1: [0, 1, 2], 2: [0, 1, 2], 3: [0, 1, 2], 4: [0]}
    >>> col_count = 5
    >>> counts = {"skip": 0.4, "s1": .33, "s2": .33}
    >>> subs = ["skip", "s1", "s2"]
    >>> conf_mat = fill_confidence_matrix(col_count, counts, subs, active, align_mat)
    >>> f"{conf_mat[1,1]:.4f}"
    '0.1100'
    >>> f"{conf_mat[2,1]:.4f}"
    '0.8800'
    >>> f"{conf_mat[1,2]:.4f}"
    '0.1980'
    >>> f"{conf_mat[2,2]:.4f}"
    '0.7920'
    >>> f"{conf_mat[1,3]:.4f}"
    '0.7920'
    >>> f"{conf_mat[2,3]:.4f}"
    '0.1980'
    """
    confidence_matrix: ConfidenceMatrix = {}

    for col_index in range(column_count):
        temp_region: List[float] = []

        for row_index in active_cells[col_index]:
            temp_region.append(align_matrix[row_index, col_index])

        temp_confidence: List[float] = confidence_cm(
            temp_region,
            subfam_counts,
            subfams,
            active_cells[col_index],
            0,
            False,
        )

        for row_index2 in range(len(active_cells[col_index])):
            confidence_matrix[
                active_cells[col_index][row_index2], col_index
            ] = temp_confidence[row_index2]

    return confidence_matrix
