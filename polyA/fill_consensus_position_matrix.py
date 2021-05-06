from typing import Dict, List, Tuple

from .matrices import ConsensusMatrixContainer
from .performance import timeit


@timeit()
def fill_consensus_position_matrix(
    row_count: int,
    column_count: int,
    start_all: int,
    subfams: List[str],
    chroms: List[str],
    starts: List[int],
    stops: List[int],
    consensus_starts: List[int],
    strands: List[str],
) -> ConsensusMatrixContainer:
    """
    Fills matrix that holds the consensus position for each subfam at that
    position in the alignment. Walks along the alignments one nucleotide at a time adding
    the consensus position to the matrix.

    At same time, fills ActiveCells.

    input:

    column_count: number of columns in alignment matrix - will be same number of
    columns in consensus_matrix
    row_count: number of rows in matrices
    start_all: min start position on chromosome/target sequences for whole alignment
    subfams: actual subfamily/consensus sequences from alignment
    chroms: actual target/chromosome sequences from alignment
    starts: start positions for all competing alignments (on target)
    stops: stop positions for all competing alignments (on target)
    consensus_starts: where alignment starts in the subfam/consensus sequence
    strands: what strand each of the alignments are on - reverse strand will count down instead of up

    output:

    ConsensusMatrixContainer

    >>> subs = ["", ".AA", "TT-"]
    >>> chrs = ["", ".AA", "TTT"]
    >>> strts = [0, 1, 0]
    >>> stps = [0, 2, 2]
    >>> con_strts = [-1, 0, 10]
    >>> strandss = ["", "+", "-"]
    >>> active, con_mat = fill_consensus_position_matrix(3, 3, 0, subs, chrs, strts, stps, con_strts, strandss)
    >>> con_mat
    {(1, 2): 0, (1, 3): 1, (2, 1): 10, (2, 2): 9, (2, 3): 9, (0, 0): 0, (0, 1): 0, (0, 2): 0}
    >>> active
    {2: [0, 1, 2], 3: [0, 1, 2], 1: [0, 2], 0: [0]}
    """
    active_cells: Dict[int, List[int]] = {}
    consensus_matrix: Dict[Tuple[int, int], int] = {}

    # start at 1 to ignore 'skip state'
    for row_index in range(1, row_count):

        if strands[row_index] == "+":
            consensus_pos = consensus_starts[row_index] - 1
            col_index: int = starts[row_index] - start_all + 1
            seq_index: int = 0

            while col_index < stops[row_index] + 1 - start_all + 1:

                # consensus pos only advances when there is not a gap in the subfam seq
                if subfams[row_index][seq_index] != "-":
                    consensus_pos += 1

                consensus_matrix[row_index, col_index] = consensus_pos

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
            seq_index2: int = 0

            while col_index2 < stops[row_index] + 1 - start_all + 1:

                if subfams[row_index][seq_index2] != "-":
                    consensus_pos2 -= 1
                consensus_matrix[row_index, col_index2] = consensus_pos2

                if chroms[row_index][seq_index2] != "-":
                    if col_index2 in active_cells:
                        active_cells[col_index2].append(row_index)
                    else:
                        active_cells[col_index2] = [0, row_index]
                    col_index2 += 1

                seq_index2 += 1

    for i in range(column_count):
        consensus_matrix[0, i] = 0
        if i not in active_cells:
            active_cells[i] = [0]

    return ConsensusMatrixContainer(active_cells, consensus_matrix)
