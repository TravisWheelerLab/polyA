import math
from typing import Dict, List, Tuple, Set

from .matrices import (
    CollapsedMatrices,
    ConsensusMatrix,
    SupportMatrix,
    SubfamAlignmentsMatrix,
)
from .performance import timeit


@timeit()
def dp_for_collapse(
    dp_rows: List[int],
    active: Dict[int, List[int]],
    support_matrix: SupportMatrix,
    non_empty: List[int],
) -> List[int]:
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
    >>> act = {0:[0,2], 2:[0,2], 3:[0,2]}
    >>> non_cols = [0, 2, 3]
    >>> sup_mat = {(0, 0): 0.5, (0, 2): 0.5, (0, 3): .1, (2, 0): 0.1, (2, 2): 0.1, (2, 3): 0.9}
    >>> p = dp_for_collapse(dp_r, act, sup_mat, non_cols)
    >>> p
    [2, 2, 2]
    """
    change: float = math.log(0.000000001)
    stay: float = math.log(1 - change)

    path: List[int] = []

    origin: Dict[Tuple[int, int], int] = {}
    dp: Dict[Tuple[int, int], float] = {}

    map_rows: Dict[int, int] = {}
    for k in range(len(dp_rows)):
        # maps active rows from matrix position to position in dp_rows
        map_rows[dp_rows[k]] = k

    # first col of prob_matrix is 0s
    for i in active[non_empty[0]]:
        dp[i, non_empty[0]] = 0.0

    # do dp and fill origin matrix
    # only touches active cells
    for i in range(1, len(non_empty)):
        curr_col = non_empty[i]
        prev_col = non_empty[i - 1]

        for row0 in active[curr_col]:
            row = map_rows[row0]

            max_value: float = -math.inf
            max_index: int = -1
            support_log: float = math.log(support_matrix[row0, curr_col])

            for row1 in active[prev_col]:
                prev_row = map_rows[row1]

                score: float = support_log + dp[row1, prev_col]

                if row == prev_row:
                    score += stay
                else:
                    score += change

                if score > max_value:
                    max_value = score
                    max_index = row1

            origin[row0, curr_col] = max_index
            dp[row0, curr_col] = max_value

    # get path from origin matrix
    # which row to start backtrace
    max_value = -math.inf
    max_row_index = 0
    for i in active[non_empty[-1]]:
        if max_value < dp[i, non_empty[-1]]:
            max_value = dp[i, non_empty[-1]]
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

    return path


@timeit()
def collapse_matrices(
    row_count: int,
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

    >>> non_cols = [1, 2, 3]
    >>> active = {1: [0, 1, 2], 2: [0, 1, 2], 3: [0, 1, 2]}
    >>> subs = ["s1", "s2", "s3"]
    >>> strandss = ["+", "-", "-"]
    >>> sup_mat = {(0, 1): 0.5, (0, 2): 0.5, (0, 3): .1, (1, 1): 0.2, (1, 2): 0.2, (1, 3): .2, (2, 1): 0.1, (2, 2): 0.1, (2, 3): 0.9}
    >>> con_mat = {(0, 1): 0, (0, 2): 1, (0, 3): 2, (1, 1): 0, (1, 2): 1, (1, 3): 2, (2, 1): 0, (2, 2): 3, (2, 3): 10}
    >>> (r, con_mat_col, strand_mat_col, sup_mat_col, sub_col, active_col, sub_col_ind, sub_aligns) = collapse_matrices(3, 1, non_cols, subs, strandss, [0,0,0], [2,2,2], active, sup_mat, con_mat)
    >>> r
    3
    >>> con_mat_col
    {(0, 1): 0, (1, 1): 0, (2, 1): 0, (0, 2): 1, (1, 2): 1, (2, 2): 3, (0, 3): 2, (1, 3): 2, (2, 3): 10}
    >>> strand_mat_col
    {(0, 1): '+', (1, 1): '-', (2, 1): '-', (0, 2): '+', (1, 2): '-', (2, 2): '-', (0, 3): '+', (1, 3): '-', (2, 3): '-'}
    >>> sup_mat_col
    {(0, 1): 0.5, (1, 1): 0.2, (2, 1): 0.1, (0, 2): 0.5, (1, 2): 0.2, (2, 2): 0.1, (0, 3): 0.1, (1, 3): 0.2, (2, 3): 0.9}
    >>> sub_col
    ['s1', 's2', 's3']
    >>> active_col
    {1: [0, 1, 2], 2: [0, 1, 2], 3: [0, 1, 2]}
    >>> sub_col_ind
    {'s1': 0, 's2': 1, 's3': 2}
    >>> sub_aligns
    {('s1', 1): (0, 0), ('s2', 1): (1, 0), ('s3', 1): (2, 0), ('s1', 2): (0, 1), ('s2', 2): (1, 1), ('s3', 2): (2, 3), ('s1', 3): (0, 2), ('s2', 3): (1, 2), ('s3', 3): (2, 10)}
    """
    # fixme - write another test for when there is overlap and mini dp is done

    consensus_matrix_collapse: Dict[Tuple[int, int], int] = {}
    strand_matrix_collapse: Dict[Tuple[int, int], str] = {}
    support_matrix_collapse: Dict[Tuple[int, int], float] = {}
    subfams_collapse: List[str] = []
    subfams_collapse_temp: Dict[
        str, int
    ] = {}  # maps subfam name to new row index
    subfams_count: Dict[str, List[int]] = {}  # original row in matrix
    active_cells_collapsed_set: Dict[int, Set[int]] = {}
    subfams_dp = set()

    subfam_alignments_collapse: SubfamAlignmentsMatrix = {}

    # assigns row num to subfams strings
    subfam_row_index: int = 0
    for i in range(row_count):
        subfam = subfams[i]
        if subfam in subfams_count:
            if subfam != "Tandem Repeat":
                subfams_count[subfam].append(i)
                # Only do the DP with non-TRs since TRs won't overlap
                subfams_dp.add(subfam)
        else:
            subfams_collapse.append(subfam)
            subfams_count[subfam] = [i]
            subfams_collapse_temp[subfam] = subfam_row_index
            subfam_row_index += 1

    for col_index in columns:
        # use set so don't add duplicates of subfams with mulitple alignments
        active_cells_collapsed_set[col_index] = set()

        for row_index in active_cells[col_index]:
            subfam = subfams[row_index]

            # if there is only one alignment for a subfam, just copy what was in the uncollapsed matrix
            # always do this for TRs
            if len(subfams_count[subfam]) == 1:
                consensus_pos = consensus_matrix[row_index, col_index]
                consensus_matrix_collapse[
                    subfams_collapse_temp[subfam], col_index
                ] = consensus_matrix[row_index, col_index]
                strand_matrix_collapse[
                    subfams_collapse_temp[subfam], col_index
                ] = strands[row_index]
                support_matrix_collapse[
                    subfams_collapse_temp[subfam], col_index
                ] = support_matrix[row_index, col_index]

                subfam_alignments_collapse[
                    subfam, col_index + start_all - 1
                ] = (row_index, consensus_pos)

            active_cells_collapsed_set[col_index].add(
                subfams_collapse_temp[subfam]
            )

    # convert all sets into lists for easier use later
    active_cells_collapse: Dict[int, List[int]] = {}
    for key in active_cells_collapsed_set:
        active_cells_collapse[key] = list(active_cells_collapsed_set[key])

    for subfam in subfams_dp:
        dp_rows = subfams_count[subfam]  # alignment rows for this subfam
        dp_active_cols: Dict[int, List[int]] = {}

        for dp_row_index in dp_rows:
            dp_col_start = starts[dp_row_index] - start_all + 1
            dp_col_stop = stops[dp_row_index] - start_all + 2
            for dp_col_index in range(dp_col_start, dp_col_stop):
                dp_active_cols.setdefault(dp_col_index, []).append(dp_row_index)

        dp_non_empty_cols = list(dp_active_cols.keys())
        dp_non_empty_cols.sort()

        dp_region_starts = []
        dp_region_stops = []
        prev_start = 0
        count = 0
        dp_region_rows = set()
        dp_active_rows: List[List[int]] = []

        # start at index 1 to look back at index 0
        for dp_col_index in range(1, len(dp_non_empty_cols)):
            dp_region_rows.update(
                dp_active_cols[dp_non_empty_cols[dp_col_index]]
            )

            if (
                dp_non_empty_cols[dp_col_index]
                != dp_non_empty_cols[dp_col_index - 1] + 1
            ):
                # new chrom region for this subfam
                # prev and cur active cols aren't sequential
                dp_region_starts.append(prev_start)
                dp_region_stops.append(dp_col_index - 1)
                prev_start = dp_col_index
                dp_active_rows.append(list(dp_region_rows))
                dp_region_rows = set()
                count += 1
        dp_region_starts.append(prev_start)
        dp_region_stops.append(len(dp_non_empty_cols) - 1)
        dp_active_rows.append(list(dp_region_rows))

        for region_index in range(len(dp_region_starts)):
            region_start = dp_region_starts[region_index]
            region_stop = dp_region_stops[region_index]
            non_empty_region = dp_non_empty_cols[region_start : region_stop + 1]
            if len(dp_active_rows[region_index]) > 1:
                # do dp
                collapse_path = dp_for_collapse(
                    dp_active_rows[region_index],
                    dp_active_cols,
                    support_matrix,
                    non_empty_region,
                )

                for dp_col_index in range(len(non_empty_region)):
                    collapse_col = non_empty_region[dp_col_index]
                    collapse_row = collapse_path[dp_col_index]
                    consensus_pos = consensus_matrix[collapse_row, collapse_col]

                    subfam_alignments_collapse[
                        subfam, collapse_col + start_all - 1
                    ] = (collapse_row, consensus_pos)
                    consensus_matrix_collapse[
                        subfams_collapse_temp[subfam], collapse_col
                    ] = consensus_pos
                    support_matrix_collapse[
                        subfams_collapse_temp[subfam], collapse_col
                    ] = support_matrix[collapse_row, collapse_col]
                    strand_matrix_collapse[
                        subfams_collapse_temp[subfam], collapse_col
                    ] = strands[collapse_row]
            else:
                # dont do dp
                for dp_col_index in range(len(non_empty_region)):
                    collapse_col = non_empty_region[dp_col_index]
                    row_index = dp_active_rows[region_index][0]

                    consensus_pos = consensus_matrix[row_index, collapse_col]
                    # only one row to pick from
                    collapse_row = dp_active_rows[region_index][0]

                    subfam_alignments_collapse[
                        subfam, collapse_col + start_all - 1
                    ] = (collapse_row, consensus_pos)

                    consensus_matrix_collapse[
                        subfams_collapse_temp[subfam], collapse_col
                    ] = consensus_pos
                    strand_matrix_collapse[
                        subfams_collapse_temp[subfam], collapse_col
                    ] = strands[row_index]
                    support_matrix_collapse[
                        subfams_collapse_temp[subfam], collapse_col
                    ] = support_matrix[row_index, collapse_col]
    row_count_update: int = len(subfams_collapse)

    return CollapsedMatrices(
        row_count_update,
        consensus_matrix_collapse,
        strand_matrix_collapse,
        support_matrix_collapse,
        subfams_collapse,
        active_cells_collapse,
        subfams_collapse_temp,
        subfam_alignments_collapse,
    )
