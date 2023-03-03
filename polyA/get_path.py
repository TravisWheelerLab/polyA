import math
from uuid import uuid4
from typing import Dict, List, Tuple

from polyA.performance import timeit


@timeit()
def get_path(
    columns: List[int],
    ids: List[str],
    subfams_collapse: List[str],
    last_column: List[float],
    active_cells_collapse: Dict[int, List[int]],
    origin_matrix: Dict[Tuple[int, int], int],
    same_subfam_change_matrix: Dict[Tuple[int, int], int],
) -> Tuple[List[int], List[str]]:
    """
    back traces through origin matrix to get most probable path through the DP matrix.
    Finds where the path switches to a different row and populates Changes and ChangesPosition.
    Reverses Changes and ChangesPosition because it's a backtrace so they are initially backwards.
    Jumps over removed/empty columns when necessary.

    assigns IDs to each col in matrices (corresponds to a nucleotide position in target/chrom
    sequence) - cols with same ID are part of same initial subfam

    input:
    temp_id: current id number being used - makes it so new ids are unique
    columns: list of non empty columns in matrices
    ids: list of ids for each column in matrices
    changes_orig: original changes from first DP trace
    changes_pos_orig: original changes_position from first DP trace
    columns_orig: original list of non empty columns before any nodes are extracted
    subfams_collapse: subfamily names for the rows
    last_column: last column in probability matrix. Used to find where to start backtrace.
    active_cells_collapse: holds which rows haves values for each column
    origin_matrix: origin matrix
    same_subfam_change_matrix: parallel to origin_matrix, if 1 - came from the same subfam, but
    got a change transition probability

    output:
    temp_id: updated current id number being used after function completes
    changes_position: which columns (positions in target/chrom seq) switch to different subfam
    changes: parallel array to changes_position - what subfam is being switches to
    updates input list ids

    >>> non_cols = [0, 1, 2, 3]
    >>> idss = ["", "", "", ""]
    >>> subs = ["s1", "s2"]
    >>> active_col = {0: [0, 1], 1: [0, 1], 2: [0, 1], 3: [0, 1]}
    >>> last_col = [-100, -10]
    >>> orig_mat = {(0, 0): 0, (1, 0): 1, (0, 1): 0, (1, 1): 0, (0, 2): 0, (1, 2): 0, (0, 3): 0, (1, 3): 1}
    >>> same_sub_mat = {}
    >>> (changes_pos, changess) = get_path(non_cols, idss, subs, last_col, active_col, orig_mat, same_sub_mat)
    >>> changes_pos
    [0, 2, 3]
    >>> changess
    ['s1', 's2']
    >>> idss[0] == idss[1]
    True
    >>> idss[2] == idss[3]
    True
    >>> idss[0] != idss[2]
    True
    """

    maxxx: float = -math.inf
    max_row_index: int = 0

    changes_position: List[int] = []
    changes: List[str] = []

    # which row to start backtrace
    for i in range(len(active_cells_collapse[columns[-1]])):
        if maxxx < last_column[i]:
            maxxx = last_column[i]
            max_row_index = active_cells_collapse[columns[-1]][i]

    prev_row_index: int = origin_matrix[max_row_index, columns[-1]]

    current_id = uuid4().hex

    ids[columns[-1]] = current_id

    # last column is a skip state pad
    changes_position.append(len(columns) - 1)

    # already added the last col, but this adds the one before cols so still start at last col
    for columns_index in range(len(columns) - 1, 1, -1):

        prev_column: int = columns[columns_index - 1]

        ids[columns[columns_index - 1]] = current_id

        if prev_row_index != origin_matrix[prev_row_index, prev_column]:
            current_id = uuid4().hex
            changes_position.append(columns_index - 1)
            changes.append(subfams_collapse[prev_row_index])
        else:
            if (prev_row_index, prev_column) in same_subfam_change_matrix:
                current_id = uuid4().hex
                changes_position.append(columns_index - 1)
                changes.append(subfams_collapse[prev_row_index])

        prev_row_index = origin_matrix[prev_row_index, prev_column]

    ids[columns[0]] = current_id
    changes_position.append(0)
    changes.append(subfams_collapse[prev_row_index])

    changes.reverse()
    changes_position.reverse()

    return changes_position, changes
