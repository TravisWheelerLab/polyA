from math import inf
from typing import Dict, List, Tuple

from polyA.calculate_score import calculate_score


# TODO: Clean up the doc string on this function
def fill_align_matrix(
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
    of nucleotides that contribute to the score - so beginning and trailing positions with less
    than chunksize nucleotides don't have lower scores

    Computes score for the first segment that does not start with a '.' by calling CalcScore()
    and from there keeps the base score and adds new chars score and subtracts old chars
    score - if a new gap is introduced, calls CalcScore() instead of adding onto base score
    ** padding of (chunksize-1)/2 added to right pad.. this way we can go all the way to the
    end of the sequence and calc alignscores without doing anything special

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
    Each cell in matrix is the alignment score of the surrounding chunksize number of nucleotides
    for that particular subfamily.

    >>> chros = ["", "..AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAT...............", "TAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAT..............."]
    >>> subs = ["", "..AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA...............", "AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA--A..............."]
    >>> strts = [0, 2, 0]
    >>> sub_mat = {"AA":1, "AT":-1, "TA":-1, "TT":1}
    >>> (c, m) = fill_align_matrix(0, 31, -5, -25, 1, subs, chros, strts, sub_mat)
    >>> c
    40
    >>> m
    {(1, 2): 31.0, (1, 3): 31.0, (1, 4): 31.0, (1, 5): 31.0, (1, 6): 31.0, (1, 7): 31.0, (1, 8): 31.0, (1, 9): 31.0, (1, 10): 31.0, (1, 11): 31.0, (1, 12): 31.0, (1, 13): 31.0, (1, 14): 31.0, (1, 15): 31.0, (1, 16): 31.0, (1, 17): 31.0, (1, 18): 31.0, (1, 19): 31.0, (1, 20): 31.0, (1, 21): 31.0, (1, 22): 31.0, (1, 23): 31.0, (1, 24): 29.0, (1, 25): 28.933333333333334, (1, 26): 28.862068965517242, (1, 27): 28.78571428571429, (1, 28): 28.703703703703702, (1, 29): 28.615384615384617, (1, 30): 28.52, (1, 31): 28.416666666666664, (1, 32): 28.304347826086953, (1, 33): 28.18181818181818, (1, 34): 28.047619047619047, (1, 35): 27.900000000000002, (1, 36): 27.736842105263158, (1, 37): 27.555555555555554, (1, 38): 27.352941176470587, (1, 39): 27.125, (2, 0): 27.125, (2, 1): 27.352941176470587, (2, 2): 27.555555555555557, (2, 3): 27.736842105263158, (2, 4): 27.9, (2, 5): 28.047619047619047, (2, 6): 28.181818181818183, (2, 7): 28.304347826086957, (2, 8): 28.416666666666668, (2, 9): 28.52, (2, 10): 28.615384615384617, (2, 11): 28.703703703703702, (2, 12): 28.785714285714285, (2, 13): 28.862068965517242, (2, 14): 28.933333333333334, (2, 15): 29.0, (2, 16): 31.0, (2, 17): 31.0, (2, 18): 31.0, (2, 19): 31.0, (2, 20): 31.0, (2, 21): 31.0, (2, 22): 5.0, (2, 23): -1.0, (2, 24): -3.0, (2, 25): -4.133333333333333, (2, 26): -5.344827586206897, (2, 27): -6.642857142857142, (2, 28): -8.037037037037036, (2, 29): -9.538461538461538, (2, 30): -11.16, (2, 31): -12.916666666666668, (2, 32): -14.82608695652174, (2, 33): -16.909090909090907, (2, 34): -19.19047619047619, (2, 35): -21.7, (2, 36): -24.473684210526315, (2, 37): -27.555555555555554, (2, 38): -31.0, (2, 39): -34.875, (0, 0): 1.0, (0, 1): 1.0, (0, 2): 1.0, (0, 3): 1.0, (0, 4): 1.0, (0, 5): 1.0, (0, 6): 1.0, (0, 7): 1.0, (0, 8): 1.0, (0, 9): 1.0, (0, 10): 1.0, (0, 11): 1.0, (0, 12): 1.0, (0, 13): 1.0, (0, 14): 1.0, (0, 15): 1.0, (0, 16): 1.0, (0, 17): 1.0, (0, 18): 1.0, (0, 19): 1.0, (0, 20): 1.0, (0, 21): 1.0, (0, 22): 1.0, (0, 23): 1.0, (0, 24): 1.0, (0, 25): 1.0, (0, 26): 1.0, (0, 27): 1.0, (0, 28): 1.0, (0, 29): 1.0, (0, 30): 1.0, (0, 31): 1.0, (0, 32): 1.0, (0, 33): 1.0, (0, 34): 1.0, (0, 35): 1.0, (0, 36): 1.0, (0, 37): 1.0, (0, 38): 1.0, (0, 39): 1.0}
    """

    num_cols: int = 0
    half_chunk: int = int((chunk_size - 1) / 2)
    align_matrix: Dict[Tuple[int, int], float] = {}

    # chunks can't start on gaps and gaps don't count when getting to the chunk_size nucls
    for i in range(1, len(chroms)):
        subfam_seq: str = subfams[i]
        chrom_seq: str = chroms[i]

        # starts at the first non '.' char, but offsets it in the matrix based on where
        # the alignments start in the seq - ex: if first alignment in the seq starts at 10,
        # will offset by 10

        seq_index: int = starts[
            i
        ] - edge_start  # place in subfam_seq and chrom_seq
        col_index = seq_index + half_chunk  # col in align_matrix
        k = half_chunk
        temp_index = seq_index
        temp_count = 0

        while temp_count < chunk_size - k:
            if chrom_seq[temp_index] != "-":
                temp_count += 1
            temp_index += 1

        offset: int = temp_index - seq_index

        # grabs first chunk - here seq_index = pos of first non '.' char
        chrom_slice: str = chrom_seq[seq_index : seq_index + offset]
        subfam_slice: str = subfam_seq[seq_index : seq_index + offset]

        # calculates score for first chunk and puts in align_matrix
        align_score: float = calculate_score(
            gap_ext, gap_init, subfam_slice, chrom_slice, "", "", sub_matrix
        )
        align_matrix[i, col_index - k] = (
            align_score * chunk_size / (chunk_size - k)
        )  # already to scale so don't need to * chunk_size and / chunk_size

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

                align_matrix[i, col_index - k] = (
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
                align_matrix[i, col_index - k] = (
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

                align_matrix[i, col_index] = (
                    align_score / num_nucls * chunk_size
                )

                if align_score == -inf:
                    del align_matrix[i, col_index]
                    break

                col_index += 1

            seq_index += 1

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
            align_matrix[i, col_index] = (
                align_score / (half_chunk + 1) * chunk_size
            )
            col_index += 1

        # max col_index is assigned to cols
        if num_cols < col_index:
            num_cols = col_index

    # assigns skip states an alignment score
    for j in range(num_cols):
        align_matrix[0, j] = float(skip_align_score)

    return num_cols, align_matrix
