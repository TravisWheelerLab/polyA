import math
from typing import Dict, List, Tuple, Optional

from .calculate_score import (
    calculate_score,
    calculate_complexity_adjusted_score,
)
from .matrices import AlignMatrix
from .performance import timeit


@timeit()
def fill_align_matrix(
    lambda_values: List[float],
    column_count: int,
    edge_start: int,
    chunk_size: int,
    gap_inits: List[float],
    gap_exts: List[float],
    skip_align_score: float,
    subfams: List[str],
    chroms: List[str],
    starts: List[int],
    stops: List[int],
    sub_matrices: List[Dict[str, int]],
    background_freqs: List[Optional[Dict[str, float]]],
) -> AlignMatrix:
    """
    Fills an alignment score matrix by calculating an alignment score (according
    to crossmatch scoring) for every segment of size `chunk_size` for all
    sequences.

    Scores are based on the surrounding `chunk_size` nucleotides in the
    alignment. Ex: column 15 in the matrix holds the aligment score for
    nucleotides at positons 0 - 30, assuming `chunk_size` = 31.

    Starting and trailing cells are treated differently. For example, column 0
    in the matrix holds the alignment score for nucleotides 0 - 15, column 1
    represents nucleotides 0 - 16, etc. Scores are weighted based on number of
    nucleotides that contribute to the score.

    This algorithm ignores all padded indices (".") in the chroms and subfams
    lists.

    Speed up by computing base score for the first segment, moves to next column
    but adding next chars score to base score and subtracting previous chars
    score from base score. Restarts and recomputes new base score when
    necessary.

    Inputs:

    everything needed for CalcScore()
    edge_start: where alignment starts on the target/chrom sequence
    chunk_size: size of nucletide chunks that are scored
    skip_align_score: alignment score to give the skip state (default = 0)
    subfams: alignments to calculate scores from
    chroms: alignments to calculate scores from
    starts: where in the target sequence the competing alignments start

    Outputs:

    align_matrix: Hash implementation of sparse 2D matrix used in pre-DP calculations.
    Key is tuple[int, int] that maps row, col to the value held in that cell of matrix. Rows
    are  subfamilies in the input alignment file, cols are nucleotide positions in the target sequence.

    >>> chros = ["", "AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAT...............", "TAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAT..............."]
    >>> subs = ["", "AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA...............", "AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA--A..............."]
    >>> strts = [0, 2, 0]
    >>> stps = [0, 39, 39]
    >>> sub_mat = [{"AA":1, "AT":-1, "TA":-1, "TT":1, "..":0}] * 3
    >>> back_freqs = [None] * 3
    >>> m = fill_align_matrix([0.0, 0.1, 0.1], 41, 0, 31, [0, -25, -25], [0, -5, -5], 1.0, subs, chros, strts, stps, sub_mat, back_freqs)
    >>> m
    {(1, 3): 3.1, (1, 4): 3.1, (1, 5): 3.1, (1, 6): 3.1, (1, 7): 3.1, (1, 8): 3.1, (1, 9): 3.1, (1, 10): 3.1, (1, 11): 3.1, (1, 12): 3.1, (1, 13): 3.1, (1, 14): 3.1, (1, 15): 3.1, (1, 16): 3.1, (1, 17): 3.1, (1, 18): 3.1, (1, 19): 3.1, (1, 20): 3.1, (1, 21): 3.1, (1, 22): 3.1, (1, 23): 3.1, (1, 24): 3.1, (1, 25): 2.9000000000000004, (1, 26): 2.8933333333333335, (1, 27): 2.8862068965517245, (1, 28): 2.878571428571429, (1, 29): 2.8703703703703702, (1, 30): 2.861538461538462, (1, 31): 2.8520000000000003, (1, 32): 2.841666666666667, (1, 33): 2.8304347826086955, (1, 34): 2.8181818181818183, (1, 35): 2.804761904761905, (1, 36): 2.7900000000000005, (1, 37): 2.7736842105263158, (1, 38): 2.7555555555555555, (1, 39): 2.735294117647059, (1, 40): 2.7125000000000004, (2, 1): 2.7125000000000004, (2, 2): 2.735294117647059, (2, 3): 2.755555555555556, (2, 4): 2.7736842105263158, (2, 5): 2.79, (2, 6): 2.804761904761905, (2, 7): 2.8181818181818183, (2, 8): 2.830434782608696, (2, 9): 2.841666666666667, (2, 10): 2.8520000000000003, (2, 11): 2.861538461538462, (2, 12): 2.8703703703703702, (2, 13): 2.8785714285714286, (2, 14): 2.8862068965517245, (2, 15): 2.8933333333333335, (2, 16): 2.9000000000000004, (2, 17): 3.1, (2, 18): 3.1, (2, 19): 3.1, (2, 20): 3.1, (2, 21): 3.1, (2, 22): 3.1, (2, 23): 0.5, (2, 24): -0.1, (2, 25): -0.30000000000000004, (2, 26): -0.41333333333333333, (2, 27): -0.5344827586206897, (2, 28): -0.6642857142857143, (2, 29): -0.8037037037037037, (2, 30): -0.9538461538461539, (2, 31): -1.116, (2, 32): -1.291666666666667, (2, 33): -1.482608695652174, (2, 34): -1.6909090909090907, (2, 35): -1.9190476190476191, (2, 36): -2.17, (2, 37): -2.447368421052632, (2, 38): -2.7555555555555555, (2, 39): -3.1, (2, 40): -3.4875000000000003, (0, 0): 1.0, (0, 1): 1.0, (0, 2): 1.0, (0, 3): 1.0, (0, 4): 1.0, (0, 5): 1.0, (0, 6): 1.0, (0, 7): 1.0, (0, 8): 1.0, (0, 9): 1.0, (0, 10): 1.0, (0, 11): 1.0, (0, 12): 1.0, (0, 13): 1.0, (0, 14): 1.0, (0, 15): 1.0, (0, 16): 1.0, (0, 17): 1.0, (0, 18): 1.0, (0, 19): 1.0, (0, 20): 1.0, (0, 21): 1.0, (0, 22): 1.0, (0, 23): 1.0, (0, 24): 1.0, (0, 25): 1.0, (0, 26): 1.0, (0, 27): 1.0, (0, 28): 1.0, (0, 29): 1.0, (0, 30): 1.0, (0, 31): 1.0, (0, 32): 1.0, (0, 33): 1.0, (0, 34): 1.0, (0, 35): 1.0, (0, 36): 1.0, (0, 37): 1.0, (0, 38): 1.0, (0, 39): 1.0, (0, 40): 1.0}
    """

    half_chunk: int = int((chunk_size - 1) / 2)
    align_matrix: Dict[Tuple[int, int], float] = {}

    # chunks can't start on gaps and gaps don't count when getting to the chunk_size nucls
    for i in range(1, len(chroms)):
        subfam_seq: str = subfams[i]
        chrom_seq: str = chroms[i]
        sub_matrix = sub_matrices[i]
        lamb = lambda_values[i]
        gap_init = gap_inits[i]
        gap_ext = gap_exts[i]

        char_complexity_adjustments = calculate_complexity_adjusted_score(
            background_freqs[i], subfam_seq, chrom_seq, lamb
        )

        # starts at the first non '.' char, but offsets it in the matrix based on where
        # the alignments start in the seq - ex: if first alignment in the seq starts at 10,
        # will offset by 10

        seq_index: int = 0  # place in subfam_seq and chrom_seq
        col_index = (
            starts[i] - edge_start + half_chunk + 1
        )  # col in align_matrix

        k = half_chunk
        temp_index = seq_index
        temp_count = 0

        while temp_count < chunk_size - k:
            # stop offset before padding starts
            if chrom_seq[temp_index] == ".":
                break
            if chrom_seq[temp_index] != "-":
                temp_count += 1
            temp_index += 1

        offset: int = temp_index - seq_index

        # normalizes for first non trailing cell
        chrom_slice: str = chrom_seq[seq_index : seq_index + offset]
        subfam_slice: str = subfam_seq[seq_index : seq_index + offset]

        # calculates score for first chunk and puts in align_matrix
        align_score: float = calculate_score(
            gap_ext,
            gap_init,
            subfam_slice,
            chrom_slice,
            "",
            "",
            sub_matrix,
            char_complexity_adjustments,
        )

        align_matrix[i, col_index - k] = lamb * (
            align_score * chunk_size / (chunk_size - k)
        )

        # scores for first part, until we get to full sized chunks
        for k in range(half_chunk - 1, -1, -1):
            if (
                chroms[i][seq_index + offset] != "-"
            ):  # if no new gap introduced, move along seq and add next nucl into score
                if subfams[i][seq_index + offset] == "-":
                    if subfams[i][seq_index + offset - 1] == "-":
                        align_score = align_score + gap_ext
                    else:
                        align_score = align_score + gap_init
                else:
                    align_score = (
                        align_score
                        + sub_matrix[
                            subfams[i][seq_index + offset]
                            + chroms[i][seq_index + offset]
                        ]
                        + char_complexity_adjustments[
                            chroms[i][seq_index + offset]
                        ]
                    )
                align_matrix[i, col_index - k] = lamb * (
                    align_score * chunk_size / (chunk_size - k)
                )

                offset += 1

            else:  # if new gap introduced, recalculate offset and call CalcScore again
                temp_index = seq_index
                temp_count = 0

                while temp_count < chunk_size - k:
                    if chrom_seq[temp_index] == ".":
                        break
                    if chrom_seq[temp_index] != "-":
                        temp_count += 1
                    temp_index += 1

                offset = temp_index - seq_index

                chrom_slice = chrom_seq[seq_index : seq_index + offset]
                subfam_slice = subfam_seq[seq_index : seq_index + offset]

                align_score = calculate_score(
                    gap_ext,
                    gap_init,
                    subfam_slice,
                    chrom_slice,
                    subfams[i][seq_index - 1],
                    chroms[i][seq_index - 1],
                    sub_matrix,
                    char_complexity_adjustments,
                )

                align_matrix[i, col_index - k] = lamb * (
                    align_score * chunk_size / (chunk_size - k)
                )

        col_index += 1
        num_nucls: int = chunk_size  # how many nucls contributed to align score

        # move to next chunk by adding next chars score and subtracting prev chars score
        while seq_index + offset < len(chrom_seq):
            temp_index = seq_index
            temp_count = 0

            # stop when you reach the last col to fill in
            if col_index > stops[i] - edge_start + 1:
                break

            # update offset if removing a gap
            if chrom_seq[seq_index] == "-":
                offset -= 1

            # skip over gap and not calc a score for the matrix
            if chrom_seq[seq_index + 1] != "-":
                # if new gap introduced, or gets rid of old gap, recalc offset, rerun CalcScore
                if (
                    chrom_seq[seq_index + offset] == "-"
                    or chrom_seq[seq_index] == "-"
                ):
                    while temp_count < chunk_size:
                        if chrom_seq[temp_index + 1] == ".":
                            break
                        if chrom_seq[temp_index + 1] != "-":
                            temp_count += 1
                        temp_index += 1

                    offset = temp_index - seq_index

                    chrom_slice = chrom_seq[
                        seq_index + 1 : seq_index + offset + 1
                    ]

                    subfam_slice = subfam_seq[
                        seq_index + 1 : seq_index + offset + 1
                    ]
                    align_score = calculate_score(
                        gap_ext,
                        gap_init,
                        subfam_slice,
                        chrom_slice,
                        subfam_seq[seq_index],
                        chrom_seq[seq_index],
                        sub_matrix,
                        char_complexity_adjustments,
                    )

                    temp_count2: int = 0
                    for nuc in chrom_slice:
                        if nuc == ".":
                            break
                        if nuc != "-":
                            temp_count2 += 1
                    num_nucls = temp_count2

                    if num_nucls <= half_chunk:
                        align_score = -math.inf

                else:
                    # align_score from previous segment - prev chars score + next chars score

                    # subtracting prev chars score
                    if subfam_seq[seq_index] == "-":
                        num_nucls -= 1
                        if subfam_seq[seq_index - 1] == "-":
                            align_score = align_score - gap_ext
                        else:
                            align_score = align_score - gap_init
                    else:
                        align_score = (
                            align_score
                            - sub_matrix[
                                subfam_seq[seq_index] + chrom_seq[seq_index]
                            ]
                            - char_complexity_adjustments[chrom_seq[seq_index]]
                        )
                        num_nucls -= 1

                    # adding next chars score
                    if subfam_seq[seq_index + offset - half_chunk] == ".":
                        break
                    elif subfam_seq[seq_index + offset] == "-":
                        num_nucls += 1
                        if subfam_seq[seq_index + offset - 1] == "-":
                            align_score = align_score + gap_ext
                        else:
                            align_score = align_score + gap_init
                    elif subfam_seq[seq_index + offset] == ".":
                        align_score = align_score
                    else:
                        align_score = (
                            align_score
                            + sub_matrix[
                                subfam_seq[seq_index + offset]
                                + chrom_seq[seq_index + offset]
                            ]
                            + char_complexity_adjustments[
                                chrom_seq[seq_index + offset]
                            ]
                        )
                        num_nucls += 1

                align_matrix[i, col_index] = lamb * (
                    align_score / num_nucls * chunk_size
                )

                if align_score == -math.inf:
                    del align_matrix[i, col_index]
                    break
                col_index += 1

            seq_index += 1

    # assigns skip states an alignment score
    # do not lambda adjust skip state score
    for j in range(column_count):
        align_matrix[0, j] = skip_align_score

    return align_matrix
