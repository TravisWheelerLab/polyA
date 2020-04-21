from math import inf
from typing import Dict, Iterable, Tuple

from calculate_chunk_score import calculate_chunk_score
from alignment import Alignment
from constants import DEFAULT_CHUNK_SIZE, DEFAULT_GAP_EXTEND, DEFAULT_GAP_START
from substitution_matrix import SubstitutionMatrix

AlignScoreMatrix = Dict[Tuple[int, int], float]
"""
TODO: Kaitlin
"""


def fill_align_score_matrix(
    alignments: Iterable[Alignment],
    gap_extend_score: int,
    gap_start_score: int,
    edge_start: int,
    substitution_matrix: SubstitutionMatrix,
    chunk_size: int = DEFAULT_CHUNK_SIZE,
) -> AlignScoreMatrix:
    """
    TODO: Kaitlin
    """
    align_score_matrix: AlignScoreMatrix = {}

    for alignment_index, alignment in enumerate(alignments):
        subfamily_sequence = alignment.subfamily_sequence
        sequence = alignment.sequence

        nuc_index = alignment.start - edge_start
        # TODO: Force chunk_size to be even?
        score_index = nuc_index + int(chunk_size / 2)

        align_score: float = 0.0

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

            align_score_matrix[alignment_index, score_index - nuc_offset] = int(
                align_score * chunk_size / (chunk_size - nuc_offset)
            )

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
                            subfamily_sequence[nuc_index + nuc_offset - 15]
                            == "."
                            or sequence[nuc_index + nuc_offset - 15] == "."
                        ):
                            align_score = -inf
                        elif subfamily_sequence[nuc_index + nuc_offset] == "-":
                            nuc_count += 1
                            if (
                                subfamily_sequence[nuc_index + nuc_offset - 1]
                                == "-"
                            ):
                                align_score += gap_extend_score
                            else:
                                align_score += gap_start_score
                        elif (
                            subfamily_sequence[nuc_index + nuc_offset] == "."
                            or sequence[nuc_index + nuc_offset] == "."
                        ):
                            align_score = align_score
                        else:
                            align_score += substitution_matrix[
                                subfamily_sequence[nuc_index + nuc_offset],
                                sequence[nuc_index + nuc_offset],
                            ]
                            nuc_count += 1

                    if align_score <= 0:
                        align_score_matrix[alignment_index, nuc_index] = 1
                    else:
                        align_score_matrix[alignment_index, nuc_index] = int(
                            align_score / nuc_count * chunk_size
                        )

                    if align_score == -inf:
                        del align_score_matrix[alignment_index, nuc_index]
                        break

            nuc_index += 1
            prev_chunk_offset = chunk_offset

        # TODO: We just need to keep track of the max nucleotide index
        if cols < nuc_index:
            cols = nuc_index

    return align_score_matrix


if __name__ == "__main__":
    import doctest

    doctest.testmod()
