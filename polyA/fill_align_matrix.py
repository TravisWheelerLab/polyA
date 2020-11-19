from math import inf
from typing import Dict, List, Tuple

from polyA.calculate_score import calculate_score


def fill_align_matrix(
    lambdaa: float,
    edge_start: int,
    chunk_size: int,
    gap_ext: int,
    gap_init: int,
    skip_align_score: float,
    subfams: List[str],
    chroms: List[str],
    starts: List[int],
    sub_matrix: Dict[str, int],
) -> Tuple[int, Dict[Tuple[int, int], float]]:
    """
    fills AlignScoreMatrix by calculating alignment score (according to crossmatch scoring)
    for every segment of size chunksize for all seqs in alignments

    Scores are of the surrounding chunksize nucleotides in the alignment. Ex: column 15 in
    matrix holds aligment score for nucleotides at positons 0 - 30. (if chunksize = 31)

    Starting and trailing cells are different - column 0 in matrix holds alignment score for
    nucleotides 0 - 15, column 1 is nucleotides 0 - 16, etc. Scores are weighted based on number
    of nucleotides that contribute to the score

    ignores all padded indices in subfams and chroms array (padded with '.')

    Speed up by computing base score for the first segment, moves to next column but adding next
    chars score to base score and subtracting previous chars score from base score. Restarts and
    recomputes new base score when necessary

    input:
    everything needed for CalcScore()
    edge_start: where alignment starts on the target/chrom sequence
    chunk_size: size of nucletide chunks that are scored
    skip_align_score: alignment score to give the skip state (default = 0)
    subfams: alignments to calculate scores from
    chroms: alignments to calculate scores from
    starts: where in the target sequence the competing alignments start

    output:
    num_cols: number of columns in the matrix
    align_matrix: Hash implementation of sparse 2D matrix used in pre-DP calculations.
    Key is tuple[int, int] that maps row, col to the value held in that cell of matrix. Rows
    are  subfamilies in the input alignment file, cols are nucleotide positions in the target sequence.

    >>> chros = ["", "..AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAT...............", "TAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAT..............."]
    >>> subs = ["", "..AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA...............", "AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA--A..............."]
    >>> strts = [0, 2, 0]
    >>> sub_mat = {"AA":1, "AT":-1, "TA":-1, "TT":1, "..":0}
    >>> (c, m) = fill_align_matrix(0.1, 0, 31, -5, -25, 1, subs, chros, strts, sub_mat)
    >>> c
    41
    >>> m
    {(1, 2): 0.7741935483870969, (1, 1): 0.6774193548387097, (1, 3): 3.1, (1, 4): 3.1, (1, 5): 3.1, (1, 6): 3.1, (1, 7): 3.1, (1, 8): 3.1, (1, 9): 3.1, (1, 10): 3.1, (1, 11): 3.1, (1, 12): 3.1, (1, 13): 3.1, (1, 14): 3.1, (1, 15): 3.1, (1, 16): 3.1, (1, 17): 3.1, (1, 18): 3.1, (1, 19): 3.1, (1, 20): 3.1, (1, 21): 3.1, (1, 22): 3.1, (1, 23): 3.1, (1, 24): 3.1, (1, 25): 2.9000000000000004, (1, 26): 2.8933333333333335, (1, 27): 2.8862068965517245, (1, 28): 2.878571428571429, (1, 29): 2.8703703703703702, (1, 30): 2.861538461538462, (1, 31): 2.8520000000000003, (1, 32): 2.841666666666667, (1, 33): 2.8304347826086955, (1, 34): 2.8181818181818183, (1, 35): 2.804761904761905, (1, 36): 2.7900000000000005, (1, 37): 2.7736842105263158, (1, 38): 2.7555555555555555, (1, 39): 2.735294117647059, (1, 40): 2.7125000000000004, (2, 1): 2.7125000000000004, (2, 2): 2.735294117647059, (2, 3): 2.755555555555556, (2, 4): 2.7736842105263158, (2, 5): 2.79, (2, 6): 2.804761904761905, (2, 7): 2.8181818181818183, (2, 8): 2.830434782608696, (2, 9): 2.841666666666667, (2, 10): 2.8520000000000003, (2, 11): 2.861538461538462, (2, 12): 2.8703703703703702, (2, 13): 2.8785714285714286, (2, 14): 2.8862068965517245, (2, 15): 2.8933333333333335, (2, 16): 2.9000000000000004, (2, 17): 3.1, (2, 18): 3.1, (2, 19): 3.1, (2, 20): 3.1, (2, 21): 3.1, (2, 22): 3.1, (2, 23): 0.5, (2, 24): -0.1, (2, 25): -0.30000000000000004, (2, 26): -0.41333333333333333, (2, 27): -0.5344827586206897, (2, 28): -0.6642857142857143, (2, 29): -0.8037037037037037, (2, 30): -0.9538461538461539, (2, 31): -1.116, (2, 32): -1.291666666666667, (2, 33): -1.482608695652174, (2, 34): -1.6909090909090907, (2, 35): -1.9190476190476191, (2, 36): -2.17, (2, 37): -2.447368421052632, (2, 38): -2.7555555555555555, (2, 39): -3.1, (2, 40): -3.4875000000000003, (0, 0): 1.0, (0, 1): 1.0, (0, 2): 1.0, (0, 3): 1.0, (0, 4): 1.0, (0, 5): 1.0, (0, 6): 1.0, (0, 7): 1.0, (0, 8): 1.0, (0, 9): 1.0, (0, 10): 1.0, (0, 11): 1.0, (0, 12): 1.0, (0, 13): 1.0, (0, 14): 1.0, (0, 15): 1.0, (0, 16): 1.0, (0, 17): 1.0, (0, 18): 1.0, (0, 19): 1.0, (0, 20): 1.0, (0, 21): 1.0, (0, 22): 1.0, (0, 23): 1.0, (0, 24): 1.0, (0, 25): 1.0, (0, 26): 1.0, (0, 27): 1.0, (0, 28): 1.0, (0, 29): 1.0, (0, 30): 1.0, (0, 31): 1.0, (0, 32): 1.0, (0, 33): 1.0, (0, 34): 1.0, (0, 35): 1.0, (0, 36): 1.0, (0, 37): 1.0, (0, 38): 1.0, (0, 39): 1.0, (0, 40): 1.0}
    """

    num_cols: int = 1
    half_chunk: int = int((chunk_size - 1) / 2)
    align_matrix: Dict[Tuple[int, int], float] = {}

    # chunks can't start on gaps and gaps don't count when getting to the chunk_size nucls
    for i in range(1, len(chroms)):
        subfam_seq: str = subfams[i]
        chrom_seq: str = chroms[i]

        # starts at the first non '.' char, but offsets it in the matrix based on where
        # the alignments start in the seq - ex: if first alignment in the seq starts at 10,
        # will offset by 10

        seq_index: int = (
            starts[i] - edge_start
        )  # place in subfam_seq and chrom_seq
        col_index = seq_index + half_chunk + 1  # col in align_matrix
        k = half_chunk
        temp_index = seq_index
        temp_count = 0

        while temp_count < chunk_size - k:
            if chrom_seq[temp_index] != "-":
                temp_count += 1
            temp_index += 1

        offset: int = temp_index - seq_index

        # add trailing cells and doesn't normalize
        for trailing in range(1, half_chunk + 1):
            if col_index - k - trailing >= 1:
                chrom_slice: str = chrom_seq[
                    seq_index : seq_index + offset - trailing
                ]
                subfam_slice: str = subfam_seq[
                    seq_index : seq_index + offset - trailing
                ]

                # calculates score for first chunk and puts in align_matrix
                align_score: float = calculate_score(
                    gap_ext,
                    gap_init,
                    subfam_slice,
                    chrom_slice,
                    "",
                    "",
                    sub_matrix,
                )

                num_nucls0 = (seq_index + offset - trailing) - seq_index + 1

                align_matrix[i, col_index - k - trailing] = lambdaa * (
                    align_score * num_nucls0 / chunk_size
                )  # already to scale so don't need to * chunk_size and / chunk_size

        # normalizes for first non trailing cell
        chrom_slice: str = chrom_seq[seq_index : seq_index + offset]
        subfam_slice: str = subfam_seq[seq_index : seq_index + offset]

        # calculates score for first chunk and puts in align_matrix
        align_score: float = calculate_score(
            gap_ext, gap_init, subfam_slice, chrom_slice, "", "", sub_matrix
        )
        align_matrix[i, col_index - k] = lambdaa * (
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
                    )

                align_matrix[i, col_index - k] = lambdaa * (
                    align_score * chunk_size / (chunk_size - k)
                )

                offset += 1

            else:  # if new gap introduced, recalculate offset and call CalcScore again
                temp_index = seq_index
                temp_count = 0

                while temp_count < chunk_size - k:
                    if chrom_seq[temp_index] != "-":
                        temp_count += 1
                    temp_index += 1

                offset = temp_index - seq_index

                chrom_slice: str = chrom_seq[seq_index : seq_index + offset]
                subfam_slice: str = subfam_seq[seq_index : seq_index + offset]

                align_score = calculate_score(
                    gap_ext,
                    gap_init,
                    subfam_slice,
                    chrom_slice,
                    subfams[i][seq_index - 1],
                    chroms[i][seq_index - 1],
                    sub_matrix,
                )
                align_matrix[i, col_index - k] = lambdaa * (
                    align_score * chunk_size / (chunk_size - k)
                )

        col_index += 1
        num_nucls: int = chunk_size  # how many nucls contributed to align score

        # move to next chunk by adding next chars score and subtracting prev chars score
        while seq_index + offset < len(chrom_seq):
            temp_index = seq_index
            temp_count = 0

            # stop when get to end of alignment and padding starts
            if chrom_seq[seq_index + 1] == ".":
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
                    )

                    temp_count2: int = 0
                    for nuc in chrom_slice:
                        if nuc == ".":
                            break
                        if nuc != "-" and nuc != ".":
                            temp_count2 += 1
                    num_nucls = temp_count2

                    if num_nucls <= half_chunk:
                        align_score = -inf

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
                        )
                        num_nucls -= 1

                    # adding next chars score
                    if subfam_seq[seq_index + offset - half_chunk] == ".":
                        align_score = -inf
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
                        )
                        num_nucls += 1

                align_matrix[i, col_index] = lambdaa * (
                    align_score / num_nucls * chunk_size
                )

                if align_score == -inf:
                    del align_matrix[i, col_index]
                    break

                col_index += 1

            seq_index += 1

        # add trailing cells and normalizes
        for trailing in range(1, half_chunk + 1):
            # if col_index - k - trailing >= 0:
            chrom_slice: str = chrom_seq[
                seq_index + trailing : seq_index + offset - 1
            ]
            subfam_slice: str = subfam_seq[
                seq_index + trailing : seq_index + offset - 1
            ]

            num_nucls2 = (seq_index + offset - 1) - (seq_index + trailing) + 1

            # calculates score for first chunk and puts in align_matrix
            align_score: float = calculate_score(
                gap_ext,
                gap_init,
                subfam_slice,
                chrom_slice,
                subfam_seq[seq_index + trailing - 1],
                chrom_seq[seq_index + trailing - 1],
                sub_matrix,
            )

            align_matrix[i, col_index - 1 + trailing] = lambdaa * (
                align_score / num_nucls2 * chunk_size
            )

        # fixes weird instance if there is a gap perfectly in the wrong place for the while loop at end
        prev_seq_index: int = seq_index
        while chrom_seq[seq_index] == "-":
            seq_index += 1

        if prev_seq_index != seq_index:
            chrom_slice = chrom_seq[seq_index::]
            subfam_slice = subfam_seq[seq_index::]

            align_score = calculate_score(
                gap_ext,
                gap_init,
                subfam_slice,
                chrom_slice,
                subfam_seq[seq_index - 1],
                chrom_seq[seq_index - 1],
                sub_matrix,
            )
            align_matrix[i, col_index] = lambdaa * (
                align_score / (half_chunk + 1) * chunk_size
            )
            col_index += 1

        # max col_index is assigned to cols
        if num_cols < col_index:
            num_cols = col_index

    # assigns skip states an alignment score
    # do not lambda adjust skip state score
    for j in range(num_cols):
        align_matrix[0, j] = skip_align_score

    # remove trailing edges that fall off end of matrix
    # can't do this during matrix construction because we don't know how many
    # cols the matrix has until the end
    for row in range(1, len(chroms)):
        for col in range(num_cols, num_cols + chunk_size + 1):
            if (row, col) in align_matrix:
                del align_matrix[row, col]

    return num_cols, align_matrix
