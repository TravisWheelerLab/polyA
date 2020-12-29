from typing import Dict, List, NamedTuple, Tuple
from math import inf, log

from polyA.matrices import CollapsedMatrices, ConsensusMatrix, SupportMatrix


def dp_for_collapse(
    subfam: str,
    dp_rows: List[int],
    support_matrix: SupportMatrix,
    columns: List[int],
    starts: List[int],
    stops: List[int],
    start_all: int,
) -> (List[int], List[int]):
    """
    When there is more than one row with the same subfam label, do a mini DP trace
    to choose which row to use in collapsed matrices

    input:
    dp_rows: rows in original matrices to be collapsed
    support_matrix: uncollapsed support matrix, use these score for dp
    columns: list of non empty cols in matrices

    output:
    path: list with labels for which row to collapse to for each column
    non_empty: list of non empty cols for this particular subfam (in the mini dp matrix)

    >>> dp_r = [0, 2]
    >>> non_cols = [0, 2, 3]
    >>> sup_mat = {(0, 0): 0.5, (0, 2): 0.5, (0, 3): .1, (1, 0): 0.2, (1, 2): 0.2, (1, 3): .2, (2, 0): 0.1, (2, 2): 0.1, (2, 3): 0.9}
    >>> (p, non) = dp_for_collapse('s', dp_r, sup_mat, non_cols, [0,0,0], [2,2,2], 1)
    >>> p
    [1, 1, 1]
    >>> non
    [0, 2, 3]
    """
    change: float = log(0.000000001)
    stay: float = log(1 - change)

    # non_empty: List[int] = []
    non_empty = set()

    path: List[int] = []

    if subfam == "Tandem Repeat":
        for c in range(len(columns)):
            col = columns[c]
            for row in dp_rows:
                if (row, col) in support_matrix and col not in non_empty:
                    non_empty.add(col)
    else:

        for i in range(len(dp_rows)):
            for j in range(
                starts[dp_rows[i]] - start_all + 1,
                stops[dp_rows[i]] - start_all + 1 + 1,
            ):
                non_empty.add(columns[j])

    non_empty = list(non_empty)

    origin: Dict[Tuple[int, int], int] = {}
    col_list: List[float] = []

    prev_col_list: List[float] = []
    # first col of prob_matrix is 0s
    for k in dp_rows:
        prev_col_list.append(0.0)

    # do dp and fill origin matrix
    for i in range(1, len(non_empty)):
        curr_col = non_empty[i]
        prev_col = non_empty[i - 1]

        col_list.clear()

        for row in range(len(dp_rows)):
            max: float = -inf
            max_index: int = -inf
            if (dp_rows[row], curr_col) in support_matrix:
                support_log: float = log(support_matrix[dp_rows[row], curr_col])

                for prev_row in range(len(prev_col_list)):

                    if (dp_rows[prev_row], prev_col) in support_matrix:
                        score: float = support_log + prev_col_list[prev_row]

                        if row == prev_row:
                            score += stay
                        else:
                            score += change

                        if score > max:
                            max = score
                            max_index = prev_row

            col_list.append(max)
            origin[row, curr_col] = max_index

        prev_col_list = col_list.copy()

    # get path from origin matrix
    # which row to start backtrace
    maxx = -inf
    max_row_index = 0
    for i in range(len(prev_col_list)):
        if maxx < prev_col_list[i]:
            maxx = prev_col_list[i]
            max_row_index = i

    prev_row_index: int = origin[max_row_index, non_empty[-1]]
    path.append(max_row_index)

    # already added the last col, but this adds the one before cols so still start at last col
    for columns_index in range(len(non_empty) - 1, 1, -1):
        prev_column: int = non_empty[columns_index - 1]
        path.append(prev_row_index)
        prev_row_index = origin[prev_row_index, prev_column]

    path.append(prev_row_index)

    path.reverse()

    return path, non_empty


def collapse_matrices(
    row_num: int,
    start_all: int,
    columns: List[int],
    subfams: List[str],
    strands: List[str],
    starts: List[int],
    stops: List[int],
    active_cells: Dict[int, List[int]],
    support_matrix: SupportMatrix,
    consensus_matrix: ConsensusMatrix,
) -> CollapsedMatrices:
    """
    In all preceding matrices, collapse and combine rows that are the same subfam.

    input:
    row_num: number of rows in matrices
    columns: list of non empty cols in matrices
    subfams: subfamily names for rows in matrices - each row is a different alignment
    strands: which strand each alignment is on
    active_cells: maps column number with rows that are active in that column
    support_matrix: uncollapsed support matrix - rows are number indices
    consensus_matrix: uncollapsed consensus matrix - rows are number indices

    output:
    CollapsedMatrices container

    >>> non_cols = [0, 2, 3]
    >>> active = {0: [0, 1, 2], 2: [0, 1, 2], 3: [0, 1, 2]}
    >>> subs = ["s1", "s2", "s1"]
    >>> strandss = ["+", "-", "-"]
    >>> sup_mat = {(0, 0): 0.5, (0, 2): 0.5, (0, 3): .1, (1, 0): 0.2, (1, 2): 0.2, (1, 3): .2, (2, 0): 0.1, (2, 2): 0.1, (2, 3): 0.9}
    >>> con_mat = {(0, 0): 0, (0, 2): 1, (0, 3): 2, (1, 0): 0, (1, 2): 1, (1, 3): 2, (2, 0): 0, (2, 2): 3, (2, 3): 10}
    >>> (r, con_mat_col, strand_mat_col, sup_mat_col, sub_col, active_col, sub_col_ind) = collapse_matrices(3, 1, non_cols, subs, strandss, [0,0,0], [2,2,2], active, sup_mat, con_mat)
    >>> r
    2
    >>> con_mat_col
    {(1, 0): 0, (1, 2): 1, (1, 3): 2, (0, 0): 0, (0, 2): 3, (0, 3): 10}
    >>> strand_mat_col
    {(1, 0): '-', (1, 2): '-', (1, 3): '-', (0, 0): '-', (0, 2): '-', (0, 3): '-'}
    >>> sup_mat_col
    {(1, 0): 0.2, (1, 2): 0.2, (1, 3): 0.2, (0, 0): 0.1, (0, 2): 0.1, (0, 3): 0.9}
    >>> sub_col
    ['s1', 's2']
    >>> active_col
    {0: [0, 1], 2: [0, 1], 3: [0, 1]}
    >>> sub_col_ind
    {'s1': 0, 's2': 1}

    """
    consensus_matrix_collapse: Dict[Tuple[int, int], int] = {}
    strand_matrix_collapse: Dict[Tuple[int, int], str] = {}
    support_matrix_collapse: Dict[Tuple[int, int], float] = {}
    subfams_collapse: List[str] = []
    subfams_collapse_temp: Dict[
        str, int
    ] = {}  # maps subfam name to new row index
    subfams_count: Dict[str, List[int]] = {}  # original row in matrix
    active_cells_collapse: Dict[int, List[int]] = {}
    subfams_dp = set()

    # assigns row num to subfams strings
    count_i: int = 0
    for i in range(row_num):
        if subfams[i] not in subfams_count:
            subfams_collapse.append(subfams[i])
            subfams_collapse_temp[subfams[i]] = count_i
            count_i += 1
            subfams_count[subfams[i]] = [i]
        else:
            subfams_count[subfams[i]].append(i)
            subfams_dp.add(subfams[i])

    for col in range(len(columns)):
        col_index: int = columns[col]
        active_cols = (
            set()
        )  # use set so don't add duplicates of subfams with mulitple alignments
        active_cells_collapse[col_index] = active_cols

        for row_index in active_cells[col_index]:
            subfam: str = subfams[row_index]

            # if there is only one alignment for a subfam, just copy what was in the uncollapsed matrix
            if len(subfams_count[subfam]) == 1:
                consensus_matrix_collapse[
                    subfams_collapse_temp[subfam], col_index
                ] = consensus_matrix[row_index, col_index]
                strand_matrix_collapse[
                    subfams_collapse_temp[subfam], col_index
                ] = strands[row_index]
                support_matrix_collapse[
                    subfams_collapse_temp[subfam], col_index
                ] = support_matrix[row_index, col_index]

            active_cells_collapse[col_index].add(subfams_collapse_temp[subfam])

    # convert all sets into lists for easier use later
    for key in active_cells_collapse:
        active_cells_collapse[key] = list(active_cells_collapse[key])

    for subfam in subfams_dp:
        collapse_path, dp_non_empty = dp_for_collapse(
            subfam,
            subfams_count[subfam],
            support_matrix,
            columns,
            starts,
            stops,
            start_all,
        )

        for i in range(len(dp_non_empty)):
            collapse_col = dp_non_empty[i]
            collapse_row = subfams_count[subfam][
                collapse_path[i]
            ]  # original row

            consensus_matrix_collapse[
                subfams_collapse_temp[subfam], collapse_col
            ] = consensus_matrix[collapse_row, collapse_col]
            support_matrix_collapse[
                subfams_collapse_temp[subfam], collapse_col
            ] = support_matrix[collapse_row, collapse_col]
            strand_matrix_collapse[
                subfams_collapse_temp[subfam], collapse_col
            ] = strands[collapse_row]

    row_num_update: int = len(subfams_collapse)

    return CollapsedMatrices(
        row_num_update,
        consensus_matrix_collapse,
        strand_matrix_collapse,
        support_matrix_collapse,
        subfams_collapse,
        active_cells_collapse,
        subfams_collapse_temp,
    )
