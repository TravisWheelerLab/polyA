from typing import Dict, List, Tuple

from polyA.matrices import ConsensusMatrixContainer


def fill_consensus_position_matrix(
    col_num: int,
    row_num: int,
    start_all: int,
    subfams: List[str],
    chroms: List[str],
    starts: List[int],
    stops: List[int],
    consensus_starts: List[int],
    strands: List[str],
) -> ConsensusMatrixContainer:
    """
    Fills parallel to AlignMatrix that holds the consensus position for each subfam at that
    position in the alignment. Walks along the alignments one nucleotide at a time adding
    the consensus position to the matrix.

    At same time, fills NonEmptyColumns and ActiveCells.

    input:
    col_num: number of columns in alignment matrix - will be same number of columns in consensus_matrix
    row_num: number of rows in matrices
    start_all: min start position on chromosome/target sequences for whole alignment
    subfams: actual subfamily/consensus sequences from alignment
    chroms: actual target/chromosome sequences from alignment
    starts: start positions for all competing alignments (on target)
    stops: stop positions for all competing alignments (on target)
    consensus_starts: where alignment starts in the subfam/consensus sequence
    strands: what strand each of the alignments are on - reverse strand will count down instead of up

    >>> subs = ["", ".AA", "TT-"]
    >>> chrs = ["", ".AA", "TTT"]
    >>> strts = [0, 1, 0]
    >>> stps = [0, 2, 2]
    >>> con_strts = [-1, 0, 10]
    >>> strandss = ["", "+", "-"]
    >>> col, active, con_mat = fill_consensus_position_matrix(3, 3, 0, subs, chrs, strts, stps, con_strts, strandss)
    >>> con_mat
    {(1, 2): 0, (1, 3): 1, (2, 1): 10, (2, 2): 9, (2, 3): 9}
    >>> col
    [0, 1, 2, 3]
    >>> active
    {2: [0, 1, 2], 3: [0, 1, 2], 1: [0, 2]}
    """
    columns = set()
    active_cells: Dict[int, List[int]] = {}
    consensus_matrix: Dict[Tuple[int, int], int] = {}

    columns.add(0)

    # start at 1 to ignore 'skip state'
    for row_index in range(1, row_num):

        if strands[row_index] == "+":
            consensus_pos = consensus_starts[row_index] - 1
            col_index: int = starts[row_index] - start_all + 1
            seq_index: int = starts[row_index] - start_all

            while col_index < stops[row_index] + 1 - start_all + 1:

                # consensus pos only advances when there is not a gap in the subfam seq
                if subfams[row_index][seq_index] != "-":
                    consensus_pos += 1

                consensus_matrix[row_index, col_index] = consensus_pos

                columns.add(col_index)

                # matrix position only advances when there is not a gap in the chrom seq
                if chroms[row_index][seq_index] != "-":
                    if col_index in active_cells:
                        active_cells[col_index].append(row_index)
                    else:
                        active_cells[col_index] = [0, row_index]
                    col_index += 1

                seq_index += 1

        else:  # reverse strand
            consensus_pos2 = consensus_starts[row_index] + 1
            col_index2: int = starts[row_index] - start_all + 1
            seq_index2: int = starts[row_index] - start_all

            while col_index2 < stops[row_index] + 1 - start_all + 1:

                if subfams[row_index][seq_index2] != "-":
                    consensus_pos2 -= 1
                consensus_matrix[row_index, col_index2] = consensus_pos2

                columns.add(col_index2)

                if chroms[row_index][seq_index2] != "-":
                    if col_index2 in active_cells:
                        active_cells[col_index2].append(row_index)
                    else:
                        active_cells[col_index2] = [0, row_index]
                    col_index2 += 1

                seq_index2 += 1

    sorted_columns: List[int] = list(columns)
    sorted_columns.sort()

    return ConsensusMatrixContainer(
        sorted_columns, active_cells, consensus_matrix
    )
