from math import inf
from typing import Dict, Iterable, Tuple

from .calculate_chunk_score import calculate_chunk_score
from .alignment import Alignment
from .constants import DEFAULT_CHUNK_SIZE, DEFAULT_GAP_EXTEND, DEFAULT_GAP_START
from .substitution_matrix import SubstitutionMatrix

AlignScoreMatrix = Dict[Tuple[int, int], int]
"""
Hash implementation of sparse 2D matrix used in pre-DP calculations. 
Key is tuple[int, int] that maps row, col to the value held in that cell of matrix. Rows 
are  subfamilies in the input alignment file, cols are nucleotide positions in the alignment.
Each cell in matrix is the alignment score of the surrounding chunksize number of nucleotides
for that particular subfamily.
"""


def fill_align_score_matrix(
    alignments: Iterable[Alignment],
    gap_extend_score: int,
    gap_start_score: int,
    edge_start: int,
    substitution_matrix: SubstitutionMatrix,
    chunk_size: int = DEFAULT_CHUNK_SIZE,
) -> Tuple[AlignScoreMatrix, int]:
    """
    
    FIXME: KAITLIN - rewrite this....
    
    fills AlignScoreMatrix by calculating alignment score (according to crossmatch scoring) 
    for every segment of size chunksize for all seqs in alignments
    
    Scores are of the surrounding chunksize nucleotides in the alignment. Ex: column 15 in 
    matrix holds aligment score for nucleotides at positons 0 - 30. 
    
    Starting and trailing cells are different - column 0 in matrix holds alignment score for
    nucleotides 0 - 15, column 1 is nucleotides 0 - 16, etc.
	
	computes score for the first segment that does not start with a '.' by calling CalcScore()
	and from there keeps the base score and adds new chars score and subtracts old chars 
	score - if a new gap is introduced, calls CalcScore() instead of adding onto base score
	
	Score are weighted based on number of nucleotides that contribute to the score - so beginning
	and trailing positions with less than chunksize nucleotides don't have lower scores
	
	TODO: KAITLIN - more robust test - check gaps on chrom seq
	
	** padding of (chunksize-1)/2 added to right pad.. this way we can go all the way to the
	end of the sequence and calc alignscores without doing anything special 
	
	>>> sub_mat = {("A", "A"): 1, ("A", "T"): -1, ("T", "A"): -1, ("T", "T"): 1}
    >>> align0 = Alignment("", 0, 2, 0, 0, 0, ["..AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAT...............", "..AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA..............."], "+")
    >>> align1 = Alignment("", 0, 0, 0, 0, 0, ["TAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAT...............", "AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA--A..............."], "+")
	>>> (matrix, col) = fill_align_score_matrix([align0, align1], -5, -25, 0, sub_mat)
	>>> matrix
	{(0, 2): 31, (0, 3): 31, (0, 4): 31, (0, 5): 31, (0, 6): 31, (0, 7): 31, (0, 8): 31, (0, 9): 31, (0, 10): 31, (0, 11): 31, (0, 12): 31, (0, 13): 31, (0, 14): 31, (0, 15): 31, (0, 16): 31, (0, 17): 31, (0, 18): 31, (0, 19): 31, (0, 20): 31, (0, 21): 31, (0, 22): 31, (0, 23): 31, (0, 24): 29, (0, 25): 28, (0, 26): 28, (0, 27): 28, (0, 28): 28, (0, 29): 28, (0, 30): 28, (0, 31): 28, (0, 32): 28, (0, 33): 28, (0, 34): 28, (0, 35): 27, (0, 36): 27, (0, 37): 27, (0, 38): 27, (0, 39): 27, (1, 0): 27, (1, 1): 27, (1, 2): 27, (1, 3): 27, (1, 4): 27, (1, 5): 28, (1, 6): 28, (1, 7): 28, (1, 8): 28, (1, 9): 28, (1, 10): 28, (1, 11): 28, (1, 12): 28, (1, 13): 28, (1, 14): 28, (1, 15): 29, (1, 16): 31, (1, 17): 31, (1, 18): 31, (1, 19): 31, (1, 20): 31, (1, 21): 31, (1, 22): 5, (1, 23): 1, (1, 24): 1, (1, 25): 1, (1, 26): 1, (1, 27): 1, (1, 28): 1, (1, 29): 1, (1, 30): 1, (1, 31): 1, (1, 32): 1, (1, 33): 1, (1, 34): 1, (1, 35): 1, (1, 36): 1, (1, 37): 1, (1, 38): 1, (1, 39): 1}
	>>> col
	40
    """
    align_score_matrix: AlignScoreMatrix = {}

    column_count = 0

    for alignment_index, alignment in enumerate(alignments):
        subfamily_sequence = alignment.subfamily_sequence
        sequence = alignment.sequence

        nuc_index = alignment.start - edge_start

        score_index = nuc_index + int((chunk_size - 1) / 2)

        align_score: float = 0.0

        chunk_offset: int = 0

        for nuc_offset in range(15, -1, -1):
            temp_index = nuc_index
            temp_count = 0

            while temp_count < chunk_size - nuc_offset:
                if sequence[temp_index] != "-":
                    temp_count += 1
                temp_index += 1

            chunk_offset = temp_index - nuc_index
            prev_chunk_offset = chunk_offset

            subfamily_slice = subfamily_sequence[
                nuc_index : nuc_index + chunk_offset
            ]
            sequence_slice = sequence[nuc_index : nuc_index + chunk_offset]

            align_score = calculate_chunk_score(
                subfamily_slice,
                sequence_slice,
                gap_extend_score,
                gap_start_score,
                substitution_matrix,
                "",
                "",
            )

            if align_score <= 0:
                align_score_matrix[
                    alignment_index, score_index - nuc_offset
                ] = 1
            else:
                align_score_matrix[
                    alignment_index, score_index - nuc_offset
                ] = int(align_score * chunk_size / (chunk_size - nuc_offset))

        score_index += 1
        nuc_count = chunk_size

        while nuc_index + chunk_offset < len(sequence):
            temp_index = nuc_index
            temp_count = 0

            while temp_count < chunk_size:
                if sequence[temp_index + 1] != "-":
                    temp_count += 1
                temp_index += 1

            chunk_offset = temp_index - nuc_index

            if sequence[nuc_index + 1] != "-":
                if (
                    sequence[nuc_index + 1] != "."
                    and subfamily_sequence[nuc_index + 1] != "."
                ):
                    if prev_chunk_offset != chunk_offset:
                        subfamily_slice = subfamily_sequence[
                            nuc_index + 1 : nuc_index + chunk_offset + 1
                        ]
                        sequence_slice = subfamily_sequence[
                            nuc_index + 1 : nuc_index + chunk_offset + 1
                        ]
                        align_score = calculate_chunk_score(
                            subfamily_slice,
                            sequence_slice,
                            gap_extend_score,
                            gap_start_score,
                            substitution_matrix,
                            "",
                            "",
                        )

                        temp_count_2 = 0
                        for nuc in sequence_slice:
                            if nuc != "-" and nuc != ".":
                                temp_count_2 += 1
                        nuc_count = temp_count_2

                        if nuc_count < 16:
                            align_score = -inf
                    else:
                        if subfamily_sequence[nuc_index] == "-":
                            nuc_count -= 1
                            if subfamily_sequence[nuc_index - 1] == "-":
                                align_score -= gap_extend_score
                            else:
                                align_score -= gap_start_score
                        else:
                            align_score -= substitution_matrix[
                                subfamily_sequence[nuc_index],
                                sequence[nuc_index],
                            ]
                            nuc_count -= 1

                        if (
                            subfamily_sequence[nuc_index + chunk_offset - 15]
                            == "."
                            or sequence[nuc_index + chunk_offset - 15] == "."
                        ):
                            align_score = -inf
                        elif (
                            subfamily_sequence[nuc_index + chunk_offset] == "-"
                        ):
                            nuc_count += 1
                            if (
                                subfamily_sequence[nuc_index + chunk_offset - 1]
                                == "-"
                            ):
                                align_score += gap_extend_score
                            else:
                                align_score += gap_start_score
                        elif (
                            subfamily_sequence[nuc_index + chunk_offset] == "."
                            or sequence[nuc_index + chunk_offset] == "."
                        ):
                            align_score = align_score
                        else:
                            align_score += substitution_matrix[
                                subfamily_sequence[nuc_index + chunk_offset],
                                sequence[nuc_index + chunk_offset],
                            ]
                            nuc_count += 1

                    if align_score <= 0:
                        align_score_matrix[alignment_index, score_index] = 1
                    else:
                        align_score_matrix[alignment_index, score_index] = int(
                            align_score / nuc_count * chunk_size
                        )

                    if align_score == -inf:
                        del align_score_matrix[alignment_index, score_index]
                        break

                score_index += 1

            nuc_index += 1
            prev_chunk_offset = chunk_offset

        #         # TODO: We just need to keep track of the max nucleotide index
        if column_count < score_index:
            column_count = score_index

    return align_score_matrix, column_count


if __name__ == "__main__":
    import doctest

    doctest.testmod()
