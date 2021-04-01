from typing import Dict, List, Tuple

from polyA.confidence_cm import confidence_cm
from polyA.matrices import ConfidenceMatrix
from polyA.performance import timeit


@timeit
def fill_confidence_matrix(
    columns: List[int],
    subfam_counts: Dict[str, float],
    subfams: List[str],
    active_cells: Dict[int, List[int]],
    align_matrix: Dict[Tuple[int, int], float],
) -> ConfidenceMatrix:
    """
    Fills confidence matrix from alignment matrix. Each column in the alignment matrix is a group of competing
    annotations that are input into confidence_cm, the output confidence values are used to populate confidence_matrix.

    input:
    everything needed for confidence_cm()
    columns: array that holds all non empty columns in align matrix
    subfam_counts: dictionary that maps subfam names to prior counts
    subfams: array of subfam names from original input alignment
    active_cells: maps col numbers to all active rows in that col
    align_matrix: matrix with alignment scores - used to calculate confidence

    output:
    confidence_matrix: Hash implementation of sparse 2D matrix used in pre-DP calculations. Key is
    tuple[int, int] that maps row, col with the value held in that cell of matrix. Rows are
    subfamilies in the input alignment file, cols are nucleotide positions in the alignment.
    Each cell in matrix is the confidence score calculated from all the alignment scores in a
    column of the AlignHash

    >>> align_mat = {(0, 0): 0, (0, 1): 100, (0, 2): 99, (1, 0): 100, (1, 1): 100, (1, 2): 100}
    >>> active = {0: [0, 1], 1: [0, 1], 2: [0, 1]}
    >>> non_cols = [0, 1, 2]
    >>> counts = {"s1": .33, "s2": .33, "s3": .33}
    >>> subs = ["s1", "s2"]
    >>> conf_mat = fill_confidence_matrix(non_cols, counts, subs, active, align_mat)
    >>> f"{conf_mat[0,0]:.4f}"
    '0.0100'
    >>> f"{conf_mat[1,0]:.4f}"
    '0.9900'
    >>> f"{conf_mat[0,1]:.4f}"
    '0.5000'
    >>> f"{conf_mat[1,1]:.4f}"
    '0.5000'
    >>> f"{conf_mat[0,2]:.4f}"
    '0.3333'
    >>> f"{conf_mat[1,2]:.4f}"
    '0.6667'
    """
    confidence_matrix: ConfidenceMatrix = {}

    for i in range(len(columns)):

        col_index: int = columns[i]
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
