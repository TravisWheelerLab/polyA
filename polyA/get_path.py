from math import inf
from typing import Dict, List, Tuple


def get_path(
    temp_id: int,
    columns: List[int],
    ids: List[int],
    changes_orig: List[str],
    changes_position_orig: List[int],
    columns_orig: List[int],
    subfams_collapse: List[str],
    last_column: List[float],
    active_cells_collapse: Dict[int, List[int]],
    origin_matrix: Dict[Tuple[int, int], int],
    same_subfam_change_matrix: Dict[Tuple[int, int], int],
) -> Tuple[int, List[int], List[str]]:
    """
    using origin matrix, back traces through the 2D array to get the subfam path (most probable
    path through the DP matrix)
    finds where the path switches to a different row and populates Changes and ChangesPosition
    reverses Changes and ChangesPosition because it's a backtrace so they are initially backwards
    jumps over removed/empty columns when necessary

    assigns IDs to each col in matrices (corresponds to a nucleotide position in target/chrom
    sequence) - cols with same ID are part of same initial subfam

    input:
    temp_id: current id number being used - makes it so new ids are unique
    columns: list of non empty columns in matrices
    ids: list of ids for each column in matrices
    changes_orig: original changes from first DP trace
    changes_pos_orig: original changes_position from first DP trace
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
    >>> idss = [0, 0, 0, 0]
    >>> subs = ["s1", "s2"]
    >>> active_col = {0: [0, 1], 1: [0, 1], 2: [0, 1], 3: [0, 1]}
    >>> last_col = [-100, -10]
    >>> orig_mat = {(0, 0): 0, (1, 0): 1, (0, 1): 0, (1, 1): 0, (0, 2): 0, (1, 2): 0, (0, 3): 0, (1, 3): 1}
    >>> same_sub_mat = {}
    >>> (temp_idd, changes_pos, changess) = get_path(1111, non_cols, idss, [], [], [], subs, last_col, active_col, orig_mat, same_sub_mat)
    >>> temp_idd
    3579
    >>> changes_pos
    [0, 2, 4]
    >>> changess
    ['s1', 's2']
    >>> idss
    [2345, 2345, 1111, 1111]
    """

    maxxx: float = -inf
    max_row_index: int = 0

    changes_position: List[int] = []
    changes: List[str] = []

    # which row to start backtrace
    for i in range(len(active_cells_collapse[columns[-1]])):
        if maxxx < last_column[i]:
            maxxx = last_column[i]
            max_row_index = active_cells_collapse[columns[-1]][i]

    prev_row_index: int = origin_matrix[max_row_index, columns[-1]]

    ids[columns[-1]] = temp_id

    #last column is a skip state pad
    changes_position.append(len(columns)-1)

    # already added the last col, but this adds the one before cols so still start at last col
    for columns_index in range(len(columns) - 1, 1, -1):

        prev_column: int = columns[columns_index - 1]
        curr_column: int = columns[columns_index]

        ids[columns[columns_index - 1]] = temp_id

        # updates the original node labels if they change when being stitched
        for i in range(len(changes_position_orig) - 1):
            if columns_orig[changes_position_orig[i]] == prev_column:
                changes_orig[i] = subfams_collapse[
                    origin_matrix[prev_row_index, curr_column]
                ]

        if prev_row_index != origin_matrix[prev_row_index, prev_column]:
            temp_id += 1234
            changes_position.append(columns_index - 1)
            changes.append(subfams_collapse[prev_row_index])
        else:
            if (prev_row_index, prev_column) in same_subfam_change_matrix:
                temp_id += 1234
                changes_position.append(columns_index - 1)
                changes.append(subfams_collapse[prev_row_index])

        prev_row_index = origin_matrix[prev_row_index, prev_column]

    ids[columns[0]] = temp_id
    changes_position.append(0)
    changes.append(subfams_collapse[prev_row_index])

    changes.reverse()
    changes_position.reverse()

    # changes ID for next round of stitching, so when starts stitching will have unique ID
    temp_id += 1234

    return temp_id, changes_position, changes
