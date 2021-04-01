from typing import Dict, List, Tuple

from .confidence_cm import confidence_cm
from .performance import timeit


@timeit
def fill_confidence_matrix_tr(
    columns: List[int],
    subfam_countss: Dict[str, float],
    subfamss: List[str],
    active_cells: Dict[int, List[int]],
    active_cells_trailing: Dict[int, List[int]],
    repeat_scores: Dict[int, float],
    align_matrix: Dict[Tuple[int, int], float],
) -> Dict[Tuple[int, int], float]:
    """
    Fills confidence matrix from alignment matrix including TR scores. Each column in the alignment matrix is a group of competing
    annotations that are input into ConfidenceCM.
    The output confidence values are used to populate confidence_matrix.

    input:
    everything needed for ConfidenceCM
    columns: array that holds all non empty columns in align matrix
    subfam_counts: dictionary that maps subfam names to prior counts
    subfams: array of subfam names from original input alignment
    active_cells: maps col numbers to all active rows in that col
    active_cells_trailing: same as above but has trailing cells on ends of alignments
    repeat_score: dict that maps col in target sequence to tandem repeat score
    align_matrix: alignment matrix - used to calculate confidence

    output:
    confidence_matrix: Hash implementation of sparse 2D matrix used in pre-DP calculations. Key is
    tuple[row, col] to value held in that cell of matrix.
    """

    confidence_matrix: Dict[Tuple[int, int], float] = {}

    # calculates confidence for alignment only columns
    # cols that feature a TR score will be replaced in next for loop
    for i in range(len(columns)):
        col_index: int = columns[i]
        temp_region: List[float] = []

        for row_index in active_cells_trailing[col_index]:
            temp_region.append(align_matrix[row_index, col_index])

        temp_confidence: List[float] = confidence_cm(
            temp_region,
            subfam_countss,
            subfamss,
            active_cells_trailing[col_index],
            0,
            False,
        )

        for row_index2 in range(len(active_cells_trailing[col_index])):
            confidence_matrix[
                active_cells_trailing[col_index][row_index2], col_index
            ] = temp_confidence[row_index2]

    # go through non empty TR columns
    for tr_col in repeat_scores:
        col_index = tr_col
        temp_region = []

        # last row in col is TR
        # assumes TRs do not overlap in a column
        for row_index in active_cells[col_index]:
            temp_region.append(align_matrix[row_index, col_index])

        temp_confidence = confidence_cm(
            temp_region,
            subfam_countss,
            subfamss,
            active_cells[col_index],
            1,
            False,
        )

        # will replace cols that had both TR and alignment scores
        for row_index2 in range(len(active_cells[col_index])):
            confidence_matrix[
                active_cells[col_index][row_index2], col_index
            ] = temp_confidence[row_index2]

    return confidence_matrix
