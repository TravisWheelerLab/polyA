from typing import Dict, List, Tuple

ConsensusPositionMatrix = Dict[Tuple[int, int], int]


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
