from typing import Dict, List, NamedTuple, Tuple

from polyA.matrices import CollapsedMatrices, ConsensusMatrix, SupportMatrix


def collapse_matrices(
    row_num: int,
    columns: List[int],
    subfams: List[str],
    strands: List[str],
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

    >>> non_cols = [0, 2, 3]
    >>> active = {0: [0, 1, 2], 2: [0, 1, 2], 3: [0, 1, 2]}
    >>> subs = ["s1", "s2", "s1"]
    >>> strandss = ["+", "-", "-"]
    >>> sup_mat = {(0, 0): 0.5, (0, 2): 0.5, (0, 3): .1, (1, 0): 0.2, (1, 2): 0.2, (1, 3): .2, (2, 0): 0.1, (2, 2): 0.1, (2, 3): 0.9}
    >>> con_mat = {(0, 0): 0, (0, 2): 1, (0, 3): 2, (1, 0): 0, (1, 2): 1, (1, 3): 2, (2, 0): 0, (2, 2): 3, (2, 3): 10}
    >>> (r, con_mat_col, strand_mat_col, sup_mat_col, sub_col, active_col, sub_col_ind) = collapse_matrices(3, non_cols, subs, strandss, active, sup_mat, con_mat)
    >>> r
    2
    >>> con_mat_col
    {(0, 0): 0, (1, 0): 0, (0, 2): 1, (1, 2): 1, (0, 3): 10, (1, 3): 2}
    >>> strand_mat_col
    {(0, 0): '+', (1, 0): '-', (0, 2): '+', (1, 2): '-', (0, 3): '-', (1, 3): '-'}
    >>> sup_mat_col
    {(0, 0): 0.5, (1, 0): 0.2, (0, 2): 0.5, (1, 2): 0.2, (0, 3): 0.9, (1, 3): 0.2}
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
    subfams_collapse_temp: Dict[str, int] = {}
    subfams_collapse: List[str] = []
    active_cells_collapse: Dict[int, List[int]] = {}

    # assigns row num to subfams strings
    count_i: int = 0
    for i in range(row_num):
        if subfams[i] not in subfams_collapse_temp:
            subfams_collapse_temp[subfams[i]] = count_i
            subfams_collapse.append(subfams[i])
            count_i += 1

    for col in range(len(columns)):
        col_index: int = columns[col]
        dup_max_support: Dict[str, float] = {}

        active_cols: List[int] = []
        active_cells_collapse[col_index] = active_cols

        # find max support score for collapsed row_nums and use that row for collapsed matrices
        for row_index in active_cells[col_index]:
            subfam: str = subfams[row_index]
            support_score: float = support_matrix[row_index, col_index]
            if (subfam) in dup_max_support:
                if support_score > dup_max_support[subfam]:
                    dup_max_support[subfam] = support_score
                    consensus_matrix_collapse[
                        subfams_collapse_temp[subfam], col_index
                    ] = consensus_matrix[row_index, col_index]
                    strand_matrix_collapse[
                        subfams_collapse_temp[subfam], col_index
                    ] = strands[row_index]
                    support_matrix_collapse[
                        subfams_collapse_temp[subfam], col_index
                    ] = support_score
            else:
                dup_max_support[subfam] = support_score
                consensus_matrix_collapse[
                    subfams_collapse_temp[subfam], col_index
                ] = consensus_matrix[row_index, col_index]
                strand_matrix_collapse[
                    subfams_collapse_temp[subfam], col_index
                ] = strands[row_index]
                support_matrix_collapse[
                    subfams_collapse_temp[subfam], col_index
                ] = support_score
                active_cells_collapse[col_index].append(
                    subfams_collapse_temp[subfam]
                )

    # update var row_nums after collapse
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
