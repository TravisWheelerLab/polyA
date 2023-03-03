import math
from typing import Dict, List, Tuple

from polyA.matrices import CollapsedMatrices
from polyA.performance import timeit


@timeit()
def fill_probability_matrix(
    same_prob_skip: float,
    same_prob: float,
    change_prob: float,
    change_prob_skip: float,
    columns: List[int],
    collapsed_matrices: CollapsedMatrices,
) -> Tuple[List[float], Dict[Tuple[int, int], int], Dict[Tuple[int, int], int]]:
    """
    Calculates the probability score matrix to find most probable path through
    the support matrix. Fills the origin matrix for easier backtrace.

    The basic algorithm is described below.

    1. Look at all i's in j-1
    2. Multiply by confidence in current cell
    3. If it comes from same i, multiply by the higher probability
    4. Else multiply by the lower probability divided by (numseqs-1) so sum of
       all probs == 1
    5. Return max

    .. note::
        All probabilities are in log space.

    Args:
        same_prob_skip: penalty given to staying in the skip state
        same_prob: penalty given to staying in the same row
        change_prob: penalty given for changing rows
        change_prob_skip: penalty given for changing rows in or out of skip state
        columns: list that holds all non empty columns in matrices
            CollapsedMatrices container

    Returns:
        Tuple:
          1. :code:`col_list`: last column of prob matrix, needed to find max to know where
             to start the backtrace.
          2. :code:`origin_matrix`: Hash implementation of sparse 2D DP matrix. This is a
             collapsed matrix. Holds which cell in previous column the probability in
             the DP matrix came from. Used when doing backtrace through the DP
             matrix.
          3. :code:`same_subfam_change_matrix`: parallel to origin_matrix, if 1 - came from
             same subfam, but got a change probability. When doing backtrace, have to
             note this is same subfam name, but different annotation.

    .. todo:: Add a larger test for this function
    """
    active_cells_collapse = collapsed_matrices.active_rows
    support_matrix_collapse = collapsed_matrices.support_matrix
    strand_matrix_collapse = collapsed_matrices.strand_matrix
    consensus_matrix_collapse = collapsed_matrices.consensus_matrix

    origin_matrix: Dict[Tuple[int, int], int] = {}
    same_subfam_change_matrix: Dict[Tuple[int, int], int] = {}
    # speed up by storing previous column (that was just calculated) in a short list for quicker access
    col_list: List[float] = []

    prev_col_list: List[float] = []
    # first col of prob_matrix is 0s
    for _ in active_cells_collapse[columns[0]]:
        prev_col_list.append(0.0)

    consensus_curr: int
    consensus_prev: int
    strand_curr: str
    strand_prev: str

    for columns_index in range(1, len(columns)):
        curr_column: int = columns[columns_index]
        prev_column: int = columns[columns_index - 1]
        col_list.clear()

        for row_index in active_cells_collapse[curr_column]:
            max_value: float = -math.inf
            max_index: int = 0
            support_log: float = math.log(
                support_matrix_collapse[row_index, curr_column]
            )
            same_subfam_change: int = 0  # if 1 - comes from the same row, but gets change prob - add to same_subfam_change_matrix and use later in GetPath()

            strand_curr = strand_matrix_collapse[
                row_index, columns[columns_index]
            ]
            consensus_curr = consensus_matrix_collapse[row_index, curr_column]

            # loop through all the rows in the previous column that have a value
            # active_cells_collapse specifies which rows have a value for each column
            temp_index: int = 0
            for prev_row_index in active_cells_collapse[prev_column]:

                score: float = support_log + prev_col_list[temp_index]
                temp_index += 1
                prob: float = change_prob

                if row_index == 0 or prev_row_index == 0:
                    prob = change_prob_skip + same_prob_skip
                    if row_index == prev_row_index:
                        prob = same_prob_skip
                else:
                    if prev_row_index == row_index:  # staying in same row
                        prob = same_prob

                        if (
                            strand_matrix_collapse[prev_row_index, prev_column]
                            != strand_curr
                        ):
                            prob = change_prob

                        elif (
                            strand_matrix_collapse[prev_row_index, prev_column]
                            == "+"
                        ):
                            if (
                                consensus_matrix_collapse[
                                    prev_row_index, prev_column
                                ]
                                > consensus_curr + 50
                            ):
                                prob = change_prob
                                same_subfam_change = 1
                        else:
                            if (
                                consensus_matrix_collapse[
                                    prev_row_index, prev_column
                                ]
                                + 50
                                < consensus_curr
                            ):
                                prob = change_prob
                                same_subfam_change = 1

                score = score + prob

                if score > max_value:
                    max_value = score
                    max_index = prev_row_index

            col_list.append(max_value)
            origin_matrix[row_index, curr_column] = max_index

            if same_subfam_change == 1 and max_index == row_index:
                same_subfam_change_matrix[row_index, curr_column] = 1

        prev_col_list = col_list.copy()

    return col_list, origin_matrix, same_subfam_change_matrix
