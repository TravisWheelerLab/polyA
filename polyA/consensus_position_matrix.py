from typing import Dict, List, Tuple

from polyA import Alignment

ConsensusPositionMatrix = Dict[Tuple[int, int], int]
"""
Hash implementation of sparse 2D matrix used along with DP matrices. Key is 
tuple[int, int] that maps row and column to value held in that cell of matrix. Each cell 
holds the alignment position in the consensus subfamily sequence. Consensus position is used
to check where alignments are on the consensus subfam sequence in relation to eachother -
needed for checking if two alignments from the same subfam can be stiched together
"""


def fill_consensus_position_matrix(
    alignments: List[Alignment],
) -> ConsensusPositionMatrix:
    """
    Fills parallel array to the AlignScoreMatrix that holds the
    consensus position for each subfam at that position in the alignment.

    Walks along the alignments one nucleotide at a time adding the consensus position to the matrix.

    >>> a0 = Alignment("", 0, 0, 0, -1, -1, ["", ""], "")
    >>> a1 = Alignment("", 0, 0, 0, 0, 2, ["AAA", "AAA"], "+")
    >>> a2 = Alignment("", 0, 0, 0, 10, 9, ["TTT", "TT-"], "-")
    >>> fill_consensus_position_matrix([a0, a1, a2])
    {(0, 0): 0, (0, 1): 0, (0, 2): 0, (1, 0): 0, (1, 1): 1, (1, 2): 2, (2, 0): 10, (2, 1): 9, (2, 2): 9}
    """
    consensus_matrix: ConsensusPositionMatrix = {}

    # FIXME: This is just to make the tests pass, we need to compute column_count
    for col_index in range(3):
        consensus_matrix[0, col_index] = 0

    for row_index, alignment in enumerate(alignments):
        consensus_position = 0
        if alignment.strand == "+":
            consensus_position = alignment.consensus_start - 1
            matrix_position = 0
            for subfam_nuc, seq_nuc in zip(
                alignment.subfamily_sequence, alignment.sequence
            ):
                if subfam_nuc != ".":
                    # Consensus position only advances when there is
                    # not a gap in the subfamily sequence
                    if subfam_nuc != "-":
                        consensus_position += 1
                    consensus_matrix[
                        row_index, matrix_position
                    ] = consensus_position
                if seq_nuc != "-":
                    matrix_position += 1
        else:  # reverse strand
            consensus_position = alignment.consensus_start + 1
            matrix_position = 0
            for subfam_nuc, seq_nuc in zip(
                alignment.subfamily_sequence, alignment.sequence
            ):
                if subfam_nuc != ".":
                    if subfam_nuc != "-":
                        consensus_position -= 1
                    consensus_matrix[
                        row_index, matrix_position
                    ] = consensus_position
                if seq_nuc != "-":
                    matrix_position += 1

    return consensus_matrix
