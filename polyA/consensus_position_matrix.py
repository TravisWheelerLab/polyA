from typing import Dict, List, Tuple

ConsensusPositionMatrix = Dict[Tuple[int, int], int]
"""
Hash implementation of sparse 2D matrix used along with DP matrices. Key is 
tuple[int, int] that maps row and column to value held in that cell of matrix. Each cell 
holds the alignment position in the consensus subfamily sequence. Consensus position is used
to check where alignments are on the consensus subfam sequence in relation to eachother -
needed for checking if two alignments from the same subfam can be stiched together
"""

def fill_consensus_position_matrix(
    subfamily_sequences: List[str],
    # TODO: Come up with a less ambiguous name for this
    sequences: List[str],
    consensus_starts: List[int],
    consensus_stops: List[int],
    column_count: int,
    row_count: int,
) -> ConsensusPositionMatrix:
    consensus_matrix: ConsensusPositionMatrix = {}
    
    """
    fills parallel array to the AlignScoreMatrix that holds the consensus position for each 
	subfam at that position in the alignment
	
	walks along the alignments one nucleotide at a time adding the consensus position to 
	the matrix
    """

    for col_index in range(column_count):
        consensus_matrix[0, col_index] = 0

    for row_index in range(row_count):
        subfam_seq = subfamily_sequences[row_index]
        seq = sequences[row_index]

        consensus_position = 0
        if strands[row_index] == "+":
            consensus_position = consensus_starts[row_index] - 1
            matrix_position = 0
            for i in range(len(subfam_seq)):
                subfam_nuc = subfam_seq[i]
                seq_nuc = seq[i]

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
        else:
            consensus_position = consensus_starts[row_index] + 1
            matrix_position = 0
            for i in range(len(subfam_seq)):
                subfam_nuc = subfam_seq[i]
                seq_nuc = seq[i]

                if subfam_nuc == ".":
                    if subfam_nuc != "-":
                        consensus_position -= 1

                    consensus_matrix[i, matrix_position] = consensus_position

                if seq_nuc != "-":
                    matrix_position += 1

    return consensus_matrix
