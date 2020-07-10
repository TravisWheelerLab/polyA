from getopt import getopt
from math import inf, log
import re
from sys import argv, stderr, stdout
from typing import Dict, List, Tuple, Union
import time
import os
import numpy as np
import json

from polyA.load_alignments import load_alignments



# -----------------------------------------------------------------------------------#
#			FUNCTIONS													   			#
# -----------------------------------------------------------------------------------#


def PrintMatrixSupport(num_col: int, start_all: int, chrom_start: int, outfile: str, matrix: Dict[Tuple[int, int], float],
                            subfams_collapse: List[str]) -> None:
    """
    prints support matrix to outfile for heatmap
    """

    with open(outfile, 'w') as out:
        out.write("\t")

        start: int = chrom_start + start_all
        j: int = 0
        while j < num_col:
            out.write(f"{start}\t")
            start += 1
            j += 1
        out.write("\n")

        for k in range(len(subfams_collapse)):
            out.write(f"{subfams_collapse[k]}\t")
            j: int = 0
            while j < num_col:
                if (k, j) in matrix:
                    out.write(str(matrix[k, j]))
                else:
                    out.write(f"{-inf}")
                out.write("\t")
                j += 1
            out.write("\n")


def PrintMatrixHash(num_col: int, num_row: int, subfams: List[str],
                    matrix: Dict[Tuple[int, int], Union[float, int, str]]) -> None:
    """
    just for debugging
    prints values inside matrices (non collapsed)
    """
    stdout.write("\t")
    j: int = 0
    while j < num_col:
        stdout.write(f"{j}\t")
        j += 1
    stdout.write("\n")

    i: int = 0
    while i < num_row:
        stdout.write(f"{subfams[i]}\t")
        j: int = 0
        while j < num_col:
            if (i, j) in matrix:
                stdout.write(f"{matrix[i, j]}")
            else:
                stdout.write(f"{-inf}")
            stdout.write("\t")
            j += 1
        stdout.write("\n")
        i += 1

def PrintMatrix(num_col: int, num_row: int, subfams: List[str],
                    matrix: np.array(Union[int, float, str])) -> None:
    """
    just for debugging
    prints values inside matrices (non collapsed)
    """
    stdout.write("\t")
    j: int = 0
    while j < num_col:
        stdout.write(f"{j}\t")
        j += 1
    stdout.write("\n")

    i: int = 0
    while i < num_row:
        stdout.write(f"{subfams[i]:<20s}\t")
        j: int = 0
        while j < num_col:
            if matrix[i * num_col + j] != 0.0:
                stdout.write(f"{matrix[i * num_col + j]}")
            else:
                stdout.write(f"{-inf:5.3g}")
            stdout.write("\t")
            j += 1
        stdout.write("\n")
        i += 1


def Edges(starts: List[int], stops: List[int]) -> Tuple[int, int]:
    """
    Find and return the min and max stop positions for the entire region
    included in the alignments.

    input:
    starts: start positions on the target sequence from the input alignment
    stops: stop positions on the target sequence from the input alignment

    output:
    minimum and maximum start and stop positions on the chromosome/target sequences for whole alignment

    >>> strt = [0, 1, 4, 7]
    >>> stp = [0, 3, 10, 9]
    >>> b, e = Edges(strt, stp)
    >>> b
    1
    >>> e
    10
    """

    min_start: int = starts[1]
    max_stop: int = stops[1]

    for i in range(1, len(starts)):
        if starts[i] < min_start:
            min_start = starts[i]
        if stops[i] > max_stop:
            max_stop = stops[i]

    return min_start, max_stop


def PadSeqs(start: List[int], stop: List[int], subfam_seqs: List[str], chrom_seqs: List[str]) -> Tuple[int, int]:
    """
    Pad out sequences with "." to allow regions where sequences do not all
    have the same start and stop positions.

    padd with an extra (chunk_size-1)/2 at the end

    input:
    start: start positions on the target sequence from the input alignment
    stop: stop positions on the target sequence from the input alignment
    subfam_seqs: actual subfamily/consensus sequences from alignment
    chrom_seqs: actual target/chromosome sequences from alignment

    output:
    updates subfam_seqs and chrom_seqs with padded sequences
    minimum and maximum start and stop positions on the chromosome/target sequences for whole alignment

    >>> strt = [0, 1, 3]
    >>> stp = [0, 1, 5]
    >>> s_seq = ['', 'a', 'aaa']
    >>> c_seq = ['', 'a', 't-t']
    >>> (b, e) = PadSeqs(strt, stp, s_seq, c_seq)
    >>> b
    1
    >>> e
    5
    >>> s_seq
    ['', 'a...................', '..aaa...............']
    >>> c_seq
    ['', 'a...................', '..t-t...............']
    """

    edge_start: int
    edge_stop: int

    (edge_start, edge_stop) = Edges(start, stop)

    for i in range(1, len(subfam_seqs)):
        left_pad: int = start[i] - edge_start
        right_pad: int = edge_stop - stop[i]

        chrom_seqs[i] = ("." * left_pad) + f"{chrom_seqs[i]}" + ("." * (right_pad + 15))
        subfam_seqs[i] = ("." * left_pad) + f"{subfam_seqs[i]}" + ("." * (right_pad + 15))

    return (edge_start, edge_stop)


def CalcScore(gap_ext: int, gap_init: int, seq1: str, seq2: str, prev_char_seq1: str,
              prev_char_seq2: str, sub_matrix: Dict[str, int]) -> float:
    """
    Calculate the score for a particular alignment between a subfamily and a target/chromsome sequence.
    Scores are calculated based on input SubstitutionMatrix, gap_ext, and gap_init.

    prev_char_seq1, prev_char_seq2 are the single nucleotides in the alignment before the chunk - if
    chunk starts with '-' these tell us to use gap_init_score or gap_extend_score as the penalty

    input:
    gap_ext: penalty to extend a gap
    gap_init: penalty to start a gap
    seq1: sequence chunk from alignment
    seq2: sequence chunk from alignment
    prev_char_seq1: character before alignment - tells if should use gap_ext or gap_init
    prev_char_seq2: character before alignment - tells if should use gap_ext or gap_init
    sub_matrix: input substitution matrix - dict that maps 2 chars being aligned to the score

    output:
    alignment score

    >>> sub_mat = {"AA":1, "AT":-1, "TA":-1, "TT":1}
    >>> CalcScore(-5, -25, "AT", "AT", "", "", sub_mat)
    2.0
    >>> CalcScore(-5, -25, "-T", "AT", "A", "A", sub_mat)
    -24.0
    >>> CalcScore(-5, -25, "-T", "AT", "-", "", sub_mat)
    -4.0
    """

    # time1: float = time.time()

    chunk_score: int = 0

    # deals with the first character of a segment being a gap character - have to look at last
    # segment to see if this is a gap init or ext
    if seq1[0] == "-":
        if prev_char_seq1 == "-":
            chunk_score += gap_ext
        else:
            chunk_score += gap_init
    elif seq2[0] == "-":
        if prev_char_seq2 == "-":
            chunk_score += gap_ext
        else:
            chunk_score += gap_init
    elif seq1[0] == ".":
        chunk_score = chunk_score
    else:
        chunk_score += sub_matrix[seq1[0] + seq2[0]]

    for j in range(1, len(seq1)):
        if seq1[j] != ".":
            if seq1[j] == "-":
                if seq1[j - 1] == "-":
                    chunk_score += gap_ext
                else:
                    chunk_score += gap_init
            elif seq2[j] == "-":
                if seq2[j - 1] == "-":
                    chunk_score += gap_ext
                else:
                    chunk_score += gap_init
            # elif seq1[j] == ".":
            #     chunk_score = chunk_score
            else:
                if (seq1[j] + seq2[j]) in sub_matrix:  # just incase the char isn't in the input substitution matrix
                    chunk_score += sub_matrix[seq1[j] + seq2[j]]
        else:
            break

    return float(chunk_score)


def FillAlignMatrix(edge_start: int, chunk_size: int, gap_ext: int, gap_init: int, skip_align_score: float,
                    subfams: List[str], chroms: List[str], starts: List[int],
                    sub_matrix: Dict[str, int]) -> Tuple[int, Dict[Tuple[int, int], float]]:
    """
    fills AlignScoreMatrix by calculating alignment score (according to crossmatch scoring)
    for every segment of size chunksize for all seqs in alignments

    Scores are of the surrounding chunksize nucleotides in the alignment. Ex: column 15 in
    matrix holds aligment score for nucleotides at positons 0 - 30. (if chunksize = 31)

    Starting and trailing cells are different - column 0 in matrix holds alignment score for
    nucleotides 0 - 15, column 1 is nucleotides 0 - 16, etc. Score are weighted based on number
    of nucleotides that contribute to the score - so beginning and trailing positions with less
    than chunksize nucleotides don't have lower scores

    computes score for the first segment that does not start with a '.' by calling CalcScore()
    and from there keeps the base score and adds new chars score and subtracts old chars
    score - if a new gap is introduced, calls CalcScore() instead of adding onto base score

    ** padding of (chunksize-1)/2 added to right pad.. this way we can go all the way to the
    end of the sequence and calc alignscores without doing anything special

    input:
    everything needed for CalcScore()
    edge_start: where alignment starts on the target/chrom sequence
    chunk_size: size of nucletide chunks that are scored

    output:
    align_matrix: Hash implementation of sparse 2D matrix used in pre-DP calculations.
    Key is tuple[int, int] that maps row, col to the value held in that cell of matrix. Rows
    are  subfamilies in the input alignment file, cols are nucleotide positions in the alignment.
    Each cell in matrix is the alignment score of the surrounding chunksize number of nucleotides
    for that particular subfamily.

    >>> chros = ["", "..AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAT...............", "TAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAT..............."]
    >>> subs = ["", "..AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA...............", "AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA--A..............."]
    >>> strts = [0, 2, 0]
    >>> sub_mat = {"AA":1, "AT":-1, "TA":-1, "TT":1}
    >>> (c, m) = FillAlignMatrix(0, 31, -5, -25, 1, subs, chros, strts, sub_mat)
    >>> c
    40
    >>> m
    {(1, 2): 31.0, (1, 3): 31.0, (1, 4): 31.0, (1, 5): 31.0, (1, 6): 31.0, (1, 7): 31.0, (1, 8): 31.0, (1, 9): 31.0, (1, 10): 31.0, (1, 11): 31.0, (1, 12): 31.0, (1, 13): 31.0, (1, 14): 31.0, (1, 15): 31.0, (1, 16): 31.0, (1, 17): 31.0, (1, 18): 31.0, (1, 19): 31.0, (1, 20): 31.0, (1, 21): 31.0, (1, 22): 31.0, (1, 23): 31.0, (1, 24): 29.0, (1, 25): 28.933333333333334, (1, 26): 28.862068965517242, (1, 27): 28.78571428571429, (1, 28): 28.703703703703702, (1, 29): 28.615384615384617, (1, 30): 28.52, (1, 31): 28.416666666666664, (1, 32): 28.304347826086953, (1, 33): 28.18181818181818, (1, 34): 28.047619047619047, (1, 35): 27.900000000000002, (1, 36): 27.736842105263158, (1, 37): 27.555555555555554, (1, 38): 27.352941176470587, (1, 39): 27.125, (2, 0): 27.125, (2, 1): 27.352941176470587, (2, 2): 27.555555555555557, (2, 3): 27.736842105263158, (2, 4): 27.9, (2, 5): 28.047619047619047, (2, 6): 28.181818181818183, (2, 7): 28.304347826086957, (2, 8): 28.416666666666668, (2, 9): 28.52, (2, 10): 28.615384615384617, (2, 11): 28.703703703703702, (2, 12): 28.785714285714285, (2, 13): 28.862068965517242, (2, 14): 28.933333333333334, (2, 15): 29.0, (2, 16): 31.0, (2, 17): 31.0, (2, 18): 31.0, (2, 19): 31.0, (2, 20): 31.0, (2, 21): 31.0, (2, 22): 5.0, (2, 23): 0.0, (2, 24): 0.0, (2, 25): 0.0, (2, 26): 0.0, (2, 27): 0.0, (2, 28): 0.0, (2, 29): 0.0, (2, 30): 0.0, (2, 31): 0.0, (2, 32): 0.0, (2, 33): 0.0, (2, 34): 0.0, (2, 35): 0.0, (2, 36): 0.0, (2, 37): 0.0, (2, 38): 0.0, (2, 39): 0.0, (0, 0): 1.0, (0, 1): 1.0, (0, 2): 1.0, (0, 3): 1.0, (0, 4): 1.0, (0, 5): 1.0, (0, 6): 1.0, (0, 7): 1.0, (0, 8): 1.0, (0, 9): 1.0, (0, 10): 1.0, (0, 11): 1.0, (0, 12): 1.0, (0, 13): 1.0, (0, 14): 1.0, (0, 15): 1.0, (0, 16): 1.0, (0, 17): 1.0, (0, 18): 1.0, (0, 19): 1.0, (0, 20): 1.0, (0, 21): 1.0, (0, 22): 1.0, (0, 23): 1.0, (0, 24): 1.0, (0, 25): 1.0, (0, 26): 1.0, (0, 27): 1.0, (0, 28): 1.0, (0, 29): 1.0, (0, 30): 1.0, (0, 31): 1.0, (0, 32): 1.0, (0, 33): 1.0, (0, 34): 1.0, (0, 35): 1.0, (0, 36): 1.0, (0, 37): 1.0, (0, 38): 1.0, (0, 39): 1.0}
    """

    time1: float = time.time()

    num_cols: int = 0

    align_matrix: Dict[Tuple[int, int], float] = {}
    # align_matrix: List[float] = []
    # row_starts: List[int] = []
    # column: List[int] = []

    # chunks can't start on gaps and gaps don't count when getting to the 30 bps
    for i in range(1, len(chroms)):
        subfam_seq: str = subfams[i]
        chrom_seq: str = chroms[i]

        # starts at the first non '.' char, but offsets it in the matrix based on where
        # the alignments start in the seq - ex: if first alignment in the seq starts at 10,
        # will offset by 10

        # calculates score for the first (chunksize-1)/2 chunks, chunk ((chunksize-1)/2)+1 is the first one that is 31nt
        seq_index: int = starts[i] - edge_start
        # col_index is the col we are in the align score matrix, $seq_index is the place in @subfam_seq and @chrom_seq
        col_index = seq_index + int((chunk_size - 1) / 2)

        k = int((chunk_size - 1) / 2)

        temp_index = seq_index
        temp_count = 0

        while temp_count < chunk_size - k:
            if chrom_seq[temp_index] != "-":
                temp_count += 1
            temp_index += 1

        offset: int = temp_index - seq_index

        # grabs first chunk - here seq_index = pos of first non '.' char
        chrom_slice: str = chrom_seq[seq_index:seq_index + offset]
        subfam_slice: str = subfam_seq[seq_index:seq_index + offset]

        # calculates score for first chunk and puts score in align_matrix
        align_score: float = CalcScore(gap_ext, gap_init, subfam_slice, chrom_slice, "", "", sub_matrix)
        align_matrix[i, col_index - k] = align_score * chunk_size / (
                    chunk_size - k)  # already to scale so don't need to * 31 and / 31
        # align_matrix.append(align_score * chunk_size / (chunk_size - k))
        # column.append(col_index - k)

        for k in range(int((chunk_size - 1) / 2) - 1, -1, -1):

            if chroms[i][seq_index + offset] != '-':

                if subfams[i][seq_index + offset] == '-':
                    if subfams[i][seq_index + offset - 1] == '-':
                        align_score = align_score + gap_ext
                    else:
                        align_score = align_score + gap_init
                else:
                    align_score = align_score + sub_matrix[
                        subfams[i][seq_index + offset] + chroms[i][seq_index + offset]]

                align_matrix[i, col_index - k] = align_score * chunk_size / (chunk_size - k)
                # align_matrix.append(align_score * chunk_size / (chunk_size - k))
                # column.append(col_index - k)

                offset += 1
            else:

                temp_index = seq_index
                temp_count = 0

                while temp_count < chunk_size - k:
                    if chrom_seq[temp_index] != "-":
                        temp_count += 1
                    temp_index += 1

                offset = temp_index - seq_index

                chrom_slice: str = chrom_seq[seq_index:seq_index + offset]
                subfam_slice: str = subfam_seq[seq_index:seq_index + offset]

                align_score = CalcScore(gap_ext, gap_init, subfam_slice, chrom_slice, subfams[i][seq_index - 1],
                                        chroms[i][seq_index - 1], sub_matrix)
                align_matrix[i, col_index - k] = align_score * chunk_size / (chunk_size - k)
                # align_matrix.append(align_score * chunk_size / (chunk_size - k))
                # column.append(col_index - k)

        col_index += 1

        num_nucls: int = chunk_size  # how many nucls contributed to align score

        # move to next chunk by adding next chars score and subtracting prev chars score
        while seq_index + offset < len(chrom_seq):
            temp_index = seq_index
            temp_count = 0

            if chrom_seq[seq_index + 1] == ".":
                break

            if chrom_seq[seq_index + 1] != "-":  #skip over gap and not calc a score for the matrix
                if chrom_seq[seq_index + offset] == '-' or chrom_seq[seq_index] == '-':  #if new gap introduced, or gets rid of old gap, recalc offset
                    while temp_count < chunk_size:
                        if chrom_seq[temp_index + 1] != "-":
                            temp_count += 1
                        temp_index += 1

                    offset = temp_index - seq_index

                    chrom_slice = chrom_seq[seq_index + 1:seq_index + offset + 1]
                    subfam_slice = subfam_seq[seq_index + 1:seq_index + offset + 1]
                    align_score = CalcScore(gap_ext, gap_init, subfam_slice, chrom_slice,
                                            subfam_seq[seq_index], chrom_seq[seq_index], sub_matrix)

                    temp_count2: int = 0
                    for nuc in chrom_slice:
                        if nuc == '.':
                            break
                        if nuc != '-' and nuc != '.':
                            temp_count2 += 1
                    num_nucls = temp_count2  # resetting num_nucls to chunk_size

                    if num_nucls <= int((chunk_size - 1) / 2):
                        align_score = -inf

                else:
                    # align_score from previous segment - prev chars score + next chars score
                    # subtracting prev  chars score - tests if its a gap in the subfam as well

                    if subfam_seq[seq_index] == "-":
                        num_nucls -= 1
                        if subfam_seq[seq_index - 1] == "-":
                            align_score = align_score - gap_ext
                        else:
                            align_score = align_score - gap_init
                    else:
                        align_score = align_score - sub_matrix[subfam_seq[seq_index] + chrom_seq[seq_index]]
                        num_nucls -= 1

                    # adding next chars score - tests if its a gap in the subfam as well
                    if subfam_seq[seq_index + offset - int((chunk_size - 1) / 2)] == ".":
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
                        align_score = align_score + sub_matrix[
                            subfam_seq[seq_index + offset] + chrom_seq[seq_index + offset]]
                        num_nucls += 1

                if align_score <= 0:
                    align_matrix[i, col_index] = 0.0
                    # align_matrix.append(0.0)
                    # column.append(col_index)
                else:
                    align_matrix[i, col_index] = align_score / num_nucls * chunk_size
                    # align_matrix.append(align_score / num_nucls * chunk_size)
                    # column.append(col_index)

                if align_score == -inf:
                    del align_matrix[i, col_index]
                    # align_matrix.pop()
                    # column.pop()
                    break

                col_index += 1

            seq_index += 1
            prev_offset = offset

        # max col_index is assigned to cols
        if num_cols < col_index:
            num_cols = col_index

        # row_starts.append(len(align_matrix))

    # assigns skip states an alignment score
    for j in range(num_cols-1, -1, -1):
        # align_matrix.insert(0,0)
        # column.insert(0,j)
        align_matrix[0, j] = float(skip_align_score)

    # for i in range(len(row_starts)):
    #     row_starts[i] += num_cols
    #
    # row_starts.insert(0, num_cols)
    # row_starts.insert(0, 0)

    print("FillAlignScoreMatrix", time.time() - time1)


    # PrintMatrixHash(num_cols, rows, Subfams, align_matrix)
    # exit()
    return (num_cols, align_matrix)


def FillConsensusPositionMatrix(col_num: int, row_num: int, start_all: int, subfams: List[str], chroms: List[str], starts: List[int], stops: List[int],
                                consensus_starts: List[int],
                                strands: List[str]) -> Dict[Tuple[int, int], int]:
    """
    Fills parallel array to the AlignScoreMatrix that holds the consensus position for each subfam
    at that position in the alignment. Walks along the alignments one nucleotide at a time adding
    the consensus position to the matrix.

    input:
    col_num: number of columns in alignment matrix - will be same number of columns in consensus_matrix
    subfams: actual subfamily/consensus sequences from alignment
    chroms: actual target/chromosome sequences from alignment
    consensus_starts: where alignment starts in the subfam/consensus sequence
    strands: what strand each of the alignments are on - reverse strand will count down instead of up

    output:
    consensus_matrix: Hash implementation of sparse 2D matrix used along with DP matrices. Key is
    tuple[int, int] that maps row and column to value help in that cell of matrix. Each cell
    holds the alignment position in the consensus subfamily sequence.

    >>> subs = ["", ".AA", "TT-"]
    >>> chrs = ["", ".AA", "TTT"]
    >>> strts = [0, 1, 0]
    >>> stps = [0, 2, 2]
    >>> con_strts = [-1, 0, 10]
    >>> strandss = ["", "+", "-"]
    >>> FillConsensusPositionMatrix(3, 3, 0, subs, chrs, strts, stps, con_strts, strandss)
    {(0, 0): 0, (0, 1): 0, (0, 2): 0, (1, 1): 0, (1, 2): 1, (2, 0): 10, (2, 1): 9, (2, 2): 9}
    """

    time1: float = time.time()

    consensus_matrix: Dict[Tuple[int, int], int] = {}

    # 0s for consensus pos of skip state
    for j in range(col_num):
        consensus_matrix[0, j] = 0

    # start at 1 to ignore 'skip state'
    for row_index in range(1, row_num):

        if strands[row_index] == "+":
            consensus_pos = consensus_starts[row_index] - 1
            col_index: int = starts[row_index] - start_all
            seq_index: int = starts[row_index] - start_all

            while col_index < stops[row_index] - start_all + 1:

                # consensus pos only advances when there is not a gap in the subfam seq
                if subfams[row_index][seq_index] != "-":
                    consensus_pos += 1

                # put consensus pos corresponding to pos in matrix in hash
                consensus_matrix[row_index, col_index] = consensus_pos

                # matrix position only advances when there is not a gap in the chrom seq
                if chroms[row_index][seq_index] != "-":
                    col_index += 1

                seq_index += 1

        else:  # reverse strand
            consensus_pos2 = consensus_starts[row_index] + 1
            col_index2: int = starts[row_index] - start_all
            seq_index2: int = starts[row_index] - start_all

            while col_index2 < stops[row_index] - start_all + 1:

                if subfams[row_index][seq_index2] != "-":
                    consensus_pos2 -= 1
                consensus_matrix[row_index, col_index2] = consensus_pos2

                if chroms[row_index][seq_index2] != "-":
                    col_index2 += 1

                seq_index2 += 1

    print("FillConsensusPositionMatrix", time.time() - time1)
    # print()

    # PrintMatrixHash(cols, rows, Subfams, consensus_matrix)
    # exit()
    return consensus_matrix


def FillColumns(num_cols: int, num_rows: int, align_matrix: Dict[Tuple[int, int], float]) -> Tuple[
    List[int], Dict[int, List[int]]]:
    """
    in input alignment some of the columns will be empty - puts all columns that are not empty into NonEmptyColumns,
    so when looping through hash I can use the vals in NonEmptyColumns - this will skip over empty columns

    input:
    num_cols: number of columns in matrices
    num_rows: number of rows in matrices
    align_matrix: alignment matrix - will look here to see which columns are empty

    output:
    columns: list of columns that are not empty

    >>> align_mat = {(0, 0): 100, (0, 2): 100, (1, 0): 100, (1, 2): 100, (2, 0): 100}
    >>> (non_cols, active) = FillColumns(3, 3, align_mat)
    >>> non_cols
    [0, 2]
    >>> active
    {0: [0, 1, 2], 2: [0, 1]}
    """

    time1: float = time.time()

    columns: List[int] = []
    active_cells: Dict[int, List[int]] = {}
    for j in range(num_cols):
        empty = 1
        i: int = 1
        active_rows: List[int] = []
        active_rows.append(0)
        for i in range(1, num_rows):
            if (i, j) in align_matrix:
                active_rows.append(i)
                empty = 0

        if not empty:
            active_cells[j] = active_rows
            columns.append(j)

    print("FillColumns", time.time() - time1)
    # exit()
    return (columns, active_cells)


def ConfidenceCM(lambdaa: float, infile: str, region: List[float], subfam_counts: Dict[str, float],
                 subfams: List[str]) -> List[float]:
    """
    computes confidence values for competing annotations using alignment scores
    Loops through the array once to find sum of 2^every_hit_score in region, then
    loops back through to calculate confidence. Converts the score to account for
    lambda before summing.

    input:
    lambdaa: lambda value for input sub_matrix (scaling factor)
    region: list of scores for competing annotations

    output:
    confidence_list: list of confidence values for competing annotations, each input alignment
    score will have one output confidence score

    >>> counts = {"s1": .33, "s2": .33, "s3": .33}
    >>> subs = ["s1", "s2", "s3"]
    >>> conf = ConfidenceCM(0.5, "infile", [2, 1, 1], counts, subs)
    >>> f"{conf[0]:.2f}"
    '0.41'
    >>> f"{conf[1]:.2f}"
    '0.29'
    >>> f"{conf[2]:.2f}"
    '0.29'
    """
    confidence_list: List[float] = []

    score_total: int = 0
    for index in range(len(region)):
        score: float = region[index]
        converted_score = score * lambdaa
        if infile:
            score_total += (2 ** converted_score) * subfam_counts[subfams[index]]
        else:
            score_total += (2 ** converted_score)

    for index in range(len(region)):
        score: float = region[index]
        converted_score = score * lambdaa
        if infile:
            confidence = ((2 ** converted_score) * subfam_counts[subfams[index]]) / score_total
        else:
            confidence = (2 ** converted_score) / score_total

        confidence_list.append(confidence)

    return confidence_list


def FillConfidenceMatrix(row_num: int, lamb: float, infilee: str, columns: List[int], subfam_countss: Dict[str, float],
                         subfamss: List[str], active_cells: Dict[int, List[int]], align_matrix: Dict[Tuple[int, int], float]) -> \
        Dict[Tuple[int, int], float]:
    """
    Fills confidence matrix from alignment matrix. Each column in the alignment matrix is a group of competing
    annotations that are input into confidence_cm, the output confidence values are used to populate the confidence
    matrix.

    input:
    row_num: number of rows in align_matrix
    lamb: lamda value for input substitution matrix
    columns: array that holds all non empty columns in align matrix
    align_matrix: alignment matrix - used to calculate confidence

    output:
    confidence_matrix: Hash implementation of sparse 2D matrix used in pre-DP calculations. Key is
    tuple[int, int] that maps row, col with the value held in that cell of matrix. Rows are
    subfamilies in the input alignment file, cols are nucleotide positions in the alignment.
    Each cell in matrix is the confidence score calculated from all the alignment scores in a
    column of the AlignHash

    >>> align_mat = {(0, 0): 0, (0, 1): 100, (0, 2): 99, (1, 0): 100, (1, 1): 100, (1, 2): 100}
    >>> active = {0: [0, 1], 1: [0, 1], 2: [0, 1]}
    >>> non_cols = [0, 1, 2]
    >>> counts = {"s1": .33, "s2": .33, "s3": .33}
    >>> subs = ["s1", "s2"]
    >>> conf_mat = FillConfidenceMatrix(2, 0.1227, "infile", non_cols, counts, subs, active, align_mat)
    >>> f"{conf_mat[0,0]:.4f}"
    '0.0002'
    >>> f"{conf_mat[1,0]:.4f}"
    '0.9998'
    >>> f"{conf_mat[0,1]:.4f}"
    '0.5000'
    >>> f"{conf_mat[1,1]:.4f}"
    '0.5000'
    >>> f"{conf_mat[0,2]:.4f}"
    '0.4788'
    >>> f"{conf_mat[1,2]:.4f}"
    '0.5212'
    """

    time1: float = time.time()

    confidence_matrix: Dict[Tuple[int, int], float] = {}

    for i in range(len(columns)):

        col_index: int = columns[i]
        temp_region: List[float] = []

        for row_index in active_cells[col_index]:
            temp_region.append(align_matrix[row_index, col_index])

        temp_confidence: List[float] = ConfidenceCM(lamb, infilee, temp_region, subfam_countss, subfamss)

        for row_index2 in range(len(active_cells[col_index])):
            confidence_matrix[active_cells[col_index][row_index2], col_index] = temp_confidence[row_index2]

    # PrintMatrixHash(cols, rows, Subfams, confidence_matrix)
    print("FillConfidenceMatrix", time.time() - time1)
    # exit()

    return confidence_matrix


def FillSupportMatrix(row_num: int, num_cols: int, chunk_size, start_all: int, columns: List[int], starts: List[int], stops:List[int], confidence_matrix: Dict[Tuple[int, int], float]) -> \
Dict[Tuple[int, int], float]:
    """
    Fills support score matrix using values in conf matrix. Score for subfam row
    at position col is sum of all confidences for subfam row for that column and
    the following chunksize-1 columns - normalized by dividing by number of segments

    score for subfam x at position i is sum of all confidences for subfam x for all segments that
    overlap position i - divided by number of segments

    input:
    row_num: number of rows in matrices
    columns: list of non empty columns in matrices
    confidence_matrix: confidence matrix - used to calculate support scores

    output:
    support_matrix: Hash implementation of sparse 2D matrix used in pre-DP calculations. Key is
    tuple[int, int] that maps row, col with the value held in that cell of matrix. Rows are
    subfamilies in the input alignment file, cols are nucleotide positions in the alignment.
    Each cell in matrix is the support score (or average confidence value) for the following
    chunksize cells in the confidence matrix.

    >>> non_cols = [0,1,2]
    >>> strts = [0, 0]
    >>> stps = [0, 2]
    >>> conf_mat = {(0, 0): 0.9, (0, 1): 0.5, (0, 2): .5, (1, 0): 0.1, (1, 1): .3, (1, 2): .1}
    >>> FillSupportMatrix(2, 3, 3, 0, non_cols, strts, stps, conf_mat)
    {(0, 0): 0.7, (0, 1): 0.6333333333333333, (0, 2): 0.5, (1, 0): 0.2, (1, 1): 0.16666666666666666, (1, 2): 0.2}
    """

    time1: float = time.time()

    support_matrix: Dict[Tuple[int, int], float] = {}

    half_chunk: int = int((chunk_size - 1) / 2)

    #skip state
    for col in range(len(columns)):
        col_index: int = columns[col]

        summ: float = 0
        num_segments: int = 0

        for sum_index in range(col_index - half_chunk, col_index + half_chunk + 1):
            if (0, sum_index) in confidence_matrix:
                num_segments += 1
                summ += confidence_matrix[0, sum_index]

        support_matrix[0, col_index] = summ / num_segments

    #rest of rows
    for row_index in range(1, row_num):

        start: int = starts[row_index] - start_all
        stop: int = stops[row_index] - start_all
        #first chunk_size/2
        left_index: int = 0
        for col_index in range(start, starts[row_index] - start_all + half_chunk):

            summ: float = 0.0
            sum_index: int = col_index - left_index
            num_segments: int = 0

            while sum_index <= col_index + half_chunk:
                summ += confidence_matrix[row_index, sum_index]
                sum_index += 1
                num_segments += 1

            support_matrix[row_index, col_index] = summ / num_segments
            left_index += 1

        # middle part where num segments is chunk_size
        #calc first chunk_size average
        col_index: int = (start + half_chunk)
        summ: float = 0.0
        sum_index: int = col_index - half_chunk
        while sum_index <= col_index + half_chunk:
            summ += confidence_matrix[row_index, sum_index]
            sum_index += 1

        support_matrix[row_index, col_index] = summ / chunk_size

        #to get next bp average subtract previous confidence, add next confidence
        for col_index in range((1 + start + half_chunk), (stop + 1) - half_chunk):
            last_index: int = col_index - half_chunk - 1
            next_index: int = col_index + half_chunk

            summ -= confidence_matrix[row_index, last_index]
            summ += confidence_matrix[row_index, next_index]

            support_matrix[row_index, col_index] = summ / chunk_size

        # last chunk_size/2
        right_index: int = half_chunk
        for col_index in range((stop + 1) - half_chunk, stops[row_index] - start_all + 1):

            summ: float = 0.0
            sum_index: int = col_index - half_chunk
            num_segments: int = 0

            while sum_index < col_index + right_index:
                summ += confidence_matrix[row_index, sum_index]
                sum_index += 1
                num_segments += 1

            support_matrix[row_index, col_index] = summ / num_segments
            right_index -= 1

    print("FillSupportMatrix", time.time() - time1)

    return support_matrix


def CollapseMatrices(row_num: int, columns: List[int], subfams: List[str], strands: List[str], active_cells: Dict[int, List[int]],
                     support_matrix: Dict[Tuple[int, int], float], consensus_matrix: Dict[Tuple[int, int], int]) -> \
        Tuple[int, Dict[Tuple[int, int], int], Dict[Tuple[int, int], str], Dict[Tuple[int, int], float], List[str],
              Dict[
                  int, List[int]], Dict[str, int]]:
    """
    collapse and combine rows that are the same subfam
    also creates a bookkeeping dict that has all the cols as keys and their values
    are arrays that hold all the active subfams in that col - used so that don't have
    to loop through all the i's just to see if a column exists

    input:
    row_num: number of rows in matrices
    columns: list of non empty cols in matrices
    subfams: subfamily names for rows in matrices - each row is a different alignment
    strands: which strand each alignment is on
    support_matrix: uncollapsed support matrix - rows are number indices
    consensus_matrix: uncollapsed consensus matrix - rows are number indices

    output:
    consensus_matrix_collapse: collapsed version
    strand_matrix_collapse: Hash implementation of sparse 2D matrix used along with DP matrices.
    Tuple[str, int] as key. String is the subfamily name of the row, int is the column in matrix.
    This is a collapsed matrix with no redundant subfamilies as rows. Each cell in matrix is
    the strand of the consensus sequence that aligned at the corresponding column position of
    the sequence.
    support_matrix_collapse: Collapsed version of Support matrix. There may be duplicate rows of the
    same subfamily, collapsing the matrices puts the duplicates into the same row. Hash
    implementation with tuple[str, int] as key. String is the subfamily name of the row, int
    is the column.
    subfams_collapse: Collapsed version of the array Subfams, any duplicate subfams are consolidated
    into one.
    active_cells_collapse: Not all rows in each column hold values. Dictionary that holds column
    number as the key, and an array of which rows hold values for that column.

    >>> non_cols = [0, 2, 3]
    >>> active = {0: [0, 1, 2], 2: [0, 1, 2], 3: [0, 1, 2]}
    >>> subs = ["s1", "s2", "s1"]
    >>> strandss = ["+", "-", "-"]
    >>> sup_mat = {(0, 0): 0.5, (0, 2): 0.5, (0, 3): .1, (1, 0): 0.2, (1, 2): 0.2, (1, 3): .2, (2, 0): 0.1, (2, 2): 0.1, (2, 3): 0.9}
    >>> con_mat = {(0, 0): 0, (0, 2): 1, (0, 3): 2, (1, 0): 0, (1, 2): 1, (1, 3): 2, (2, 0): 0, (2, 2): 3, (2, 3): 10}
    >>> (r, con_mat_col, strand_mat_col, sup_mat_col, sub_col, active_col, sub_col_ind) = CollapseMatrices(3, non_cols, subs, strandss, active, sup_mat, con_mat)
    >>> r
    2
    >>> con_mat_col
    {(0, 0): 0, (1, 0): 0, (0, 2): 1, (1, 2): 1, (0, 3): 10, (1, 3): 2}
    >>> strand_mat_col
    {(0, 0): '+', (1, 0): '-', (0, 2): '+', (1, 2): '-', (0, 3): '-', (1, 3): '-'}
    >>> sup_mat_col
    {(0, 0): 0.5, (1, 0): 0.2, (0, 2): 0.5, (1, 2): 0.2, (0, 3): 0.9, (1, 3): 0.2}
    >>> sub_col
    ['s1', 's2']
    >>> active_col
    {0: [0, 1], 2: [0, 1], 3: [0, 1]}
    >>> sub_col_ind
    {'s1': 0, 's2': 1}

    """

    time1: float = time.time()

    consensus_matrix_collapse: Dict[Tuple[int, int], int] = {}
    strand_matrix_collapse: Dict[Tuple[int, int], str] = {}
    support_matrix_collapse: Dict[Tuple[int, int], float] = {}
    subfams_collapse_temp: Dict[str, int] = {}
    subfams_collapse: List[str] = []
    active_cells_collapse: Dict[int, List[int]] = {}

    #assigns row num to subfams strings
    count_i: int = 0
    for i in range(row_num):
        if subfams[i] not in subfams_collapse_temp:
            subfams_collapse_temp[subfams[i]] = count_i
            subfams_collapse.append(subfams[i])
            count_i += 1

    for col in range(len(columns)):
        col_index: int = columns[col]
        dup_max_support: Dict[str, float] = {}

        active_cols: List[int] = []
        active_cells_collapse[col_index] = active_cols

        # find max support score for collapsed row_nums and use that row for collapsed matrices
        # for row_index in range(row_num):
        for row_index in active_cells[col_index]:
            if (subfams[row_index]) in dup_max_support:
                if support_matrix[row_index, col_index] > dup_max_support[subfams[row_index]]:
                    dup_max_support[subfams[row_index]] = support_matrix[row_index, col_index]
                    consensus_matrix_collapse[subfams_collapse_temp[subfams[row_index]], col_index] = consensus_matrix[
                            row_index, col_index]
                    strand_matrix_collapse[subfams_collapse_temp[subfams[row_index]], col_index] = strands[row_index]
                    support_matrix_collapse[subfams_collapse_temp[subfams[row_index]], col_index] = support_matrix[row_index, col_index]
            else:
                dup_max_support[subfams[row_index]] = support_matrix[row_index, col_index]
                consensus_matrix_collapse[subfams_collapse_temp[subfams[row_index]], col_index] = consensus_matrix[row_index, col_index]
                strand_matrix_collapse[subfams_collapse_temp[subfams[row_index]], col_index] = strands[row_index]
                support_matrix_collapse[subfams_collapse_temp[subfams[row_index]], col_index] = support_matrix[row_index, col_index]
                active_cells_collapse[col_index].append(subfams_collapse_temp[subfams[row_index]])

    # update var row_nums after collapse
    row_num_update: int = len(subfams_collapse)

    print("CollapseMatrices", time.time() - time1)

    return row_num_update, consensus_matrix_collapse, strand_matrix_collapse, support_matrix_collapse, subfams_collapse, active_cells_collapse, subfams_collapse_temp

#FIXME - can we make this faster? numpy arrays instead of hashes?
def FillProbabilityMatrix(same_prob_skip: float, same_prob: float, change_prob: float, change_prob_skip: float,
                          columns: List[int], subfams_collapse: List[str],
                          active_cells_collapse: Dict[int, List[int]],
                          support_matrix_collapse: Dict[Tuple[int, int], float],
                          strand_matrix_collapse: Dict[Tuple[int, int], str],
                          consensus_matrix_collapse: Dict[Tuple[int, int], int]) -> Tuple[
    Dict[Tuple[int, int], float], Dict[Tuple[int, int], int], Dict[Tuple[int, int], int]]:
    """
    Fills in the probability score matrix from the support matrix. Also fills
    the origin matrix for convenience.

    We omit filling the first column because we want its probabilities to be zero anyway.

    The basic algorithm is described below. All calculations happen in log space.
    look at all i's in j-1
        mult by confidence in current cell
        if comes from same i, mult by higher prob
        else - mult by lower prob /(numseqs-1) -> so sum of all probs == 1
     return max

     NOTE:
        all probabilities are in log space
        row indices are subfam names

    input:
    same_prob_skip: penalty given to staying in the skip state
    same_prob: penalty given to staying in the same row
    change_prob: penalty given for channging rows
    change_prob_skip: penalty given for changinr rows in or out of skip state
    columns: list that holds all non empty columns in matrices
    active_cells_collapse: holds which rows have values for all columns - using this making it so don't have to
    loop through all cells in the previous column when filling a cell, just loop though cells that hold a value
    support_matrix_collapse: probabilites are calculated form support scores
    strand_matrix_collapse: used when testing if some collapsed seqs are continuous - can't be if they aren't
    on same strand
    consensus_matrix_collapse: used when testing if some collapsed seqs are continuous - positions in consensus
    sequence have to be contiguous

    output:
    probability_matrix: Hash implementation of sparse 2D DP matrix. This is a collapsed matrix. Holds DP
    probabilies when finding most probable path through the matrix. Dict that hols subfam name as row and col number
    and maps it to probability in cell
    origin_matrix: Hash implementation of sparse 2D DP matrix. This is a collapsed matrix. Holds which cell in
    previous column the probability in the DP matrix came from. Used when doing backtrace through the DP matrix.
    Dict that hols subfam name as row and col number and maps it to value in cell

    TODO: larger test needed for this function
    """

    time1: float = time.time()

    # time_total: float = 0.0

    prob_matrix: Dict[Tuple[int, int], float] = {}
    origin_matrix: Dict[Tuple[int, int], int] = {}
    same_subfam_change_matrix: Dict[Tuple[int, int], int] = {}

    prev_col_list = []
    # fill first col of prob_matrix with 0s
    for k in active_cells_collapse[columns[0]]:
        prob_matrix[k, columns[0]] = 0
        prev_col_list.append(0.0)

    # print(active_cells_collapse[columns[0]])
    for columns_index in range(1, len(columns)):
        curr_column: int = columns[columns_index]
        prev_column: int = columns[columns_index - 1]

        #speed up by storing previous column (that was just calculated) in a short list for quicker access
        col_list = []

        for row_index in active_cells_collapse[curr_column]:
            max: float = -inf
            max_index: int = 0
            support_log: float = log(support_matrix_collapse[row_index, curr_column])
            same_subfam_change: int = 0  # if 1 - comes from the same row, but gets change prob - add to same_subfam_change_matrix and use later in GetPath()
            # row = prev_row_index
            # loop through all the rows in the previous column that have a value - active_cells_collapse
            # specifies which rows have a value for each column

            temp_index: int = 0
            for prev_row_index in active_cells_collapse[prev_column]:
                # time_s = time.time()
                #just this lookup in prob_matrix is 40% of the runtime of the whole function
                score: float = support_log + prev_col_list[temp_index]#prob_matrix[prev_row_index, prev_column]
                temp_index += 1
                # time_total += time.time() - time_s
                prob: float = 0.0
                # time_s = time.time()
                #this chunk is 35% of the runtime
                if prev_row_index == row_index:  # staying in same row
                    prob = same_prob
                    if prev_row_index == 0:  # staying in skip
                        prob = same_prob_skip

                    # because rows are collapsed, if the consensus seqs are NOT contiguous - treat as if they are not the same row and get the jump penalty
                    if strand_matrix_collapse[prev_row_index, prev_column] != strand_matrix_collapse[
                        row_index, columns[columns_index]]:
                        # if not on same strand give change prob
                        prob = change_prob
                    else:  # if on same strand
                        if strand_matrix_collapse[prev_row_index, prev_column] == '+':
                            if consensus_matrix_collapse[prev_row_index, prev_column] > \
                                    consensus_matrix_collapse[row_index, curr_column] + 50:
                                prob = change_prob
                                same_subfam_change = 1
                        if strand_matrix_collapse[prev_row_index, prev_column] == '-':
                            if consensus_matrix_collapse[prev_row_index, prev_column] + 50 < \
                                    consensus_matrix_collapse[row_index, curr_column]:
                                prob = change_prob
                                same_subfam_change = 1

                else:  # jumping rows
                    prob = change_prob
                    if prev_row_index == 0 or row_index == 0:  # jumping in or out of skip
                        prob = change_prob_skip

                # time_total += time.time() - time_s

                score = score + prob

                if score > max:
                    max = score
                    max_index = prev_row_index

            col_list.append(max)
            prob_matrix[row_index, curr_column] = max
            origin_matrix[row_index, curr_column] = max_index

            if same_subfam_change == 1 and max_index == row_index:
                same_subfam_change_matrix[row_index, curr_column] = 1

        prev_col_list = col_list.copy()

    print("FillProbabilityMatrix", time.time() - time1)

    # exit()
    return (prob_matrix, origin_matrix, same_subfam_change_matrix)


def GetPath(num_col: int, temp_id: int, columns: List[int], ids: List[int], changes_orig: List[str],
            changes_position_orig: List[int], columns_orig: List[int], subfams_collapse: List[str],
            active_cells_collapse: Dict[int, List[int]], prob_matrix: Dict[Tuple[int, int], float],
            origin_matrix: Dict[Tuple[int, int], int], same_subfam_change_matrix: Dict[Tuple[int, int], int]) -> Tuple[
    int, List[int], List[str]]:
    """
    using origin matrix, back traces through the 2D array to get the subfam path (most probable
    path through the DP matrix)
    finds where the path switches to a different row and populates Changes and ChangesPosition
    reverses Changes and ChangesPosition because it's a backtrace so they are initially backwards
    jumps over removed/empty columns when necessary

    assigns IDs to each col in matrices (corresponds to a nucleotide position in target/chrom
    sequence) - cols with same ID are part of same initial subfam

    input:
    num_col: number of cols in matrices
    temp_id: current id number being used - makes it so new ids are unique
    columns: list of non empty columns in matrices
    ids: list of ids for each column in matrices
    subfams: subfamily names for the rows
    active_cells_collapse: holds which rows haves values for each column
    prob_matrix: probability matrix
    origin_matrix: origin matrix

    output:
    temp_id: updated current id number being used after function completes
    changes_position: which columns (positions in target/chrom seq) switch to different subfam
    changes: parallel array to changes_position - what subfam is being switches to
    updates input list ids

    >>> non_cols = [0, 1, 2, 3]
    >>> idss = [0, 0, 0, 0]
    >>> subs = ["s1", "s2"]
    >>> active_col = {0: [0, 1], 1: [0, 1], 2: [0, 1], 3: [0, 1]}
    >>> prob_mat = {(0, 0): 0, (1, 0): 0, (0, 1): 0, (1, 1): 0, (0, 2): 1, (1, 2): 1, (0, 3): -100, (1, 3): -10}
    >>> orig_mat = {(0, 0): 0, (1, 0): 1, (0, 1): 0, (1, 1): 0, (0, 2): 0, (1, 2): 0, (0, 3): 0, (1, 3): 1}
    >>> same_sub_mat = {}
    >>> (temp_idd, changes_pos, changess) = GetPath(4, 1111, non_cols, idss, [], [], [], subs, active_col, prob_mat, orig_mat, same_sub_mat)
    >>> temp_idd
    3579
    >>> changes_pos
    [0, 2, 4]
    >>> changess
    ['s1', 's2']
    >>> idss
    [2345, 2345, 1111, 1111]
    """

    time1: float = time.time()

    maxxx: float = -inf
    max_row_index: int = 0

    changes_position: List[int] = []
    changes: List[str] = []

    for i in active_cells_collapse[num_col - 1]:
        if maxxx < prob_matrix[i, num_col - 1]:
            maxxx = prob_matrix[i, num_col - 1]
            max_row_index = i

    prev_row_index: int = origin_matrix[max_row_index, num_col - 1]

    ids[columns[-1]] = temp_id

    changes_position.append(len(columns))

    # already added the last col, but this adds the one before cols so still start at last col
    for columns_index in range(len(columns) - 1, 1, -1):

        prev_column: int = columns[columns_index - 1]
        curr_column: int = columns[columns_index]

        ids[columns[columns_index - 1]] = temp_id

        # updates the original node labels if they change when being stitched
        for i in range(len(changes_position_orig) - 1):
            if columns_orig[changes_position_orig[i]] == prev_column:
                changes_orig[i] = subfams_collapse[origin_matrix[prev_row_index, curr_column]]

        if prev_row_index != origin_matrix[prev_row_index, prev_column]:
            temp_id += 1234
            changes_position.append(columns_index - 1)
            changes.append(subfams_collapse[prev_row_index])
        else:
            if (prev_row_index, prev_column) in same_subfam_change_matrix:
                temp_id += 1234
                changes_position.append(columns_index - 1)
                changes.append(subfams_collapse[prev_row_index])

        prev_row_index = origin_matrix[prev_row_index, prev_column]

    ids[columns[0]] = temp_id
    changes_position.append(0)
    changes.append(subfams_collapse[prev_row_index])

    changes.reverse()
    changes_position.reverse()

    # changes ID for next round of stitching, so when starts stitching will have unique ID
    temp_id += 1234

    print("GetPath", time.time() - time1)

    return (temp_id, changes_position, changes)


def PrintChanges(columns: List[int], changes: List[str], changes_position: List[int]) -> None:
    """
    just for debugging
    prints out the changes positions and subfams in the matrices, use for each iteration of the node extraction
    """
    i: int = 0
    while i < len(changes):
        stdout.write(str(columns[changes_position[i]]))
        stdout.write("\t")
        stdout.write(f"{changes[i]}\n")
        i = i + 1


def PrintNodeConfidence(nodes: int, changes: List[str], subfams_collapse: Dict[str, int],
                        node_confidence: Dict[Tuple[str, int], float]) -> None:
    """
    used for debugging
    prints out matrix that holds node confidence
    """
    for i in range(nodes):
        stdout.write(f"{changes[i]} ")
    stdout.write("\n")

    for subfam in subfams_collapse:
        stdout.write(f"{subfam} ")
        for j in range(nodes):
            if (subfam, j) in node_confidence:
                stdout.write(f"{node_confidence[subfam, j]} ")
            else:
                stdout.write(f"-inf ")
        stdout.write("\n")


def FillNodeConfidence(nodes: int, start_all: int, gap_init: int, gap_ext: int, lamb: float, infilee: str, columns: List[int],
                       starts: List[int], stops: List[int], changes_position: List[int], subfams: List[str], subfam_seqs: List[str], chrom_seqs: List[str],
                       subfam_countss: Dict[str, float], sub_matrix: Dict[str, int]) -> Dict[Tuple[str, int], float]:
    """
    finds completing annoations, alignment scores and confidence values for each node
    identified in GetPath()

    first fills matrix with node alignment scores, then reuses matrix for confidence scores

    input:
    all input needed for CalcScore() and ConfidenceCM()
    nodes: number of nodes
    changes_pos: node boundaries

    output:
    node_confidence: Hash implementation of sparse 2D matrix that holds confidence values
    for whole nodes. Used during stitching process. Tuple[str, int] is key that maps a subfamily
    and node number to a confidence score.

    >>> non_cols = [0, 1, 2, 3, 4, 5, 6, 7, 8, 9]
    >>> strts = [0, 0, 0]
    >>> stps = [10, 10, 10]
    >>> change_pos = [0, 3, 7, 10]
    >>> names = ["skip", "n1", "n2"]
    >>> s_seqs = ['', 'AAA-TTTTT-', 'TTTTTTTTTT']
    >>> c_seqs = ['', 'TTTTTTTTTT', 'TTTTTTTTTT']
    >>> counts = {"skip": .33, "n1": .33, "n2": .33}
    >>> sub_mat = {"AA":1, "AT":-1, "TA":-1, "TT":1}
    >>> node_conf = FillNodeConfidence(3, 0, -25, -5, 0.1227, "infile", non_cols, strts, stps, change_pos, names, s_seqs, c_seqs, counts, sub_mat)
    >>> node_conf
    {('skip', 0): 0.0, ('n1', 0): 0.3751243838973974, ('n2', 0): 0.6248756161026026, ('skip', 1): 0.0, ('n1', 1): 0.09874227070127324, ('n2', 1): 0.9012577292987267, ('skip', 2): 0.0, ('n1', 2): 0.09874227070127327, ('n2', 2): 0.9012577292987267}
    """

    time1: float = time.time()

    node_confidence_temp: List[float] = [0.0 for _ in range(len(subfams) * nodes)]
    node_confidence: Dict[Tuple[str, int], float] = {}

    #holds all chrom seqs and the offset needed to get to the correct position in the chrom
    #seq - remember gaps in chrom seq are skipped over for matrix position
    chrom_seq_offset: Dict[int, int] = {}

    #first node
    begin_node0: int = columns[changes_position[0]]

    for subfam_index0 in range(1, len(subfams)):

        count: int = 0
        for i in range(begin_node0, columns[changes_position[1]]-columns[changes_position[0]]):
            if chrom_seqs[subfam_index0][i] == '-':
                count += 1

        chrom_seq_offset[subfam_index0] = count

        end_node0: int = columns[changes_position[1]] + count

        subfam0: str = subfam_seqs[subfam_index0][begin_node0:end_node0]
        chrom0: str = chrom_seqs[subfam_index0][begin_node0:end_node0]

        align_score0: float = 0.0
        #if whole alignment is padding - don't run CalcScore
        if end_node0 >= starts[subfam_index0]-start_all and begin_node0 <= stops[subfam_index0] - start_all:
            align_score0 = CalcScore(gap_ext, gap_init, subfam0, chrom0, '', '', sub_matrix)
        node_confidence_temp[subfam_index0 * nodes + 0] = align_score0

    #middle nodes
    for node_index in range(1, nodes-1):

        for subfam_index in range(1, len(subfams)):

            begin_node: int = columns[changes_position[node_index]] + chrom_seq_offset[subfam_index]

            count: int = 0
            for i in range(begin_node, begin_node + columns[changes_position[node_index + 1]] - columns[changes_position[node_index]]):
                if chrom_seqs[subfam_index][i] == '-':
                    count += 1

            chrom_seq_offset[subfam_index] = chrom_seq_offset[subfam_index] + count

            end_node: int = columns[changes_position[node_index+1]] + chrom_seq_offset[subfam_index] + count

            lastprev_subfam: str = subfam_seqs[subfam_index][begin_node-1]
            lastprev_chrom: str = chrom_seqs[subfam_index][begin_node - 1]
            subfam: str = subfam_seqs[subfam_index][begin_node:end_node]
            chrom: str = chrom_seqs[subfam_index][begin_node:end_node]

            align_score: float = 0.0
            # if whole alignment is padding - don't run CalcScore
            if end_node >= starts[subfam_index] - start_all and begin_node <= stops[subfam_index] - start_all:
                align_score = CalcScore(gap_ext, gap_init, subfam, chrom, lastprev_subfam, lastprev_chrom, sub_matrix)

            node_confidence_temp[subfam_index * nodes + node_index] = align_score


    # does last node
    for subfam_index2 in range(1, len(subfams)):
        begin_node2: int = columns[changes_position[-2]] + chrom_seq_offset[subfam_index2]

        count: int = 0
        for i in range(begin_node2,
                       begin_node2 + columns[changes_position[-1] - 1] - columns[changes_position[-2]]):
            if chrom_seqs[subfam_index2][i] == '-':
                count += 1

        chrom_seq_offset[subfam_index2] = chrom_seq_offset[subfam_index2] + count

        end_node2: int = columns[changes_position[-1] - 1] + chrom_seq_offset[subfam_index2] + count

        lastprev_subfam2: str = subfam_seqs[subfam_index2][begin_node2 - 1]
        lastprev_chrom2: str = chrom_seqs[subfam_index2][begin_node2 - 1]

        subfam2: str = subfam_seqs[subfam_index2][begin_node2:end_node2+1]
        chrom2: str = chrom_seqs[subfam_index2][begin_node2:end_node2+1]

        align_score2: float = 0.0
        # if whole alignment is padding - don't run CalcScore
        if end_node2 >= starts[subfam_index2] - start_all and begin_node0 <= stops[subfam_index2] - start_all:
            align_score2 = CalcScore(gap_ext, gap_init, subfam2, chrom2, lastprev_subfam2, lastprev_chrom2,
                                                sub_matrix)

        node_confidence_temp[subfam_index2 * nodes + nodes - 1] = align_score2


    # reuse same matrix and compute confidence scores for the nodes
    for node_index4 in range(nodes):
        temp: List[float] = []
        for row_index in range(1, len(subfams)):
            temp.append(node_confidence_temp[row_index * nodes + node_index4])

        confidence_temp: List[float] = ConfidenceCM(lamb, infilee, temp, subfam_countss, subfams)

        for row_index2 in range(len(confidence_temp)):
            node_confidence_temp[(row_index2+1) * nodes + node_index4] = confidence_temp[row_index2]

    # collapse node_confidence down same way supportmatrix is collapsed - all seqs of
    # the same subfam are put in the same row
    # not a sparse hash - holds the 0s, but I think this is okay because it won't ever
    # be a very large matrix, and this way we don't have to test if anything exists

    for node_index5 in range(nodes):
        for row_index3 in range(len(subfams)):
            if (subfams[row_index3], node_index5) in node_confidence:
                node_confidence[subfams[row_index3], node_index5] += node_confidence_temp[
                    row_index3 * nodes + node_index5]
            else:
                node_confidence[subfams[row_index3], node_index5] = node_confidence_temp[
                    row_index3 * nodes + node_index5]

    print("FillNodeConfidence", time.time() - time1)

    return node_confidence


def PrintPathGraph(nodes: int, changes: List[str], path_graph: List[int]) -> None:
    """
    used for debugging
    prints out the matrix that holds the path graph - 1 means there is an edge, 0 means no edge
    """
    stdout.write(" ")
    for i in range(nodes):
        stdout.write(f"{changes[i]} ")
    stdout.write("\n")

    for i in range(nodes):
        for j in range(nodes):
            stdout.write(f"{path_graph[i * nodes + j]}\t")
        stdout.write(f"{changes[i]}\n")
    stdout.write("\n")


def FillPathGraph(nodes: int, columns: List[int], changes: List[str], changes_position: List[int], consensus_matrix_collapse: Dict[Tuple[int, int], int],
                  strand_matrix_collapse: Dict[Tuple[int, int], str], node_confidence: Dict[Tuple[str, int], float], subfams_collapse_index: Dict[str, int]) -> \
        List[int]:
    """
    finds alternative paths through the nodes - used for stitching to find nested elements
    and stitch back together elements that have been inserted into

    input:
    nodes: number of nodes
    columns: list with all non empty columns in matrices
    changes: list of node boundaries
    changes_position: list of subfams identified for each node
    subfams_collapse: list of all non duplicate subfams
    consensus_matrix_collapse: matrix holding alignment position in consensus/subfam seqs -
    used to identify whether nodes are spacially close enough for an alternative edge
    strand_matrix_collapse: matrix holding strand for alignments - nodes can only be connected
    with alternative edge if on same strand
    node_confidence: competing annoation confidence for nodes - nodes can only be connected
    with alternative edge if subfam identified for sink node is a competing annotation in
    source node and is above a certain confidence

    output:
    path_graph: 2D matrix. Graph used during stitching. Maps nodes to all other nodes and holds
    values for if there is an alternative edge between the nodes.

    >>> non_cols = [0, 1, 2, 3, 4, 5, 6, 7, 8, 9]
    >>> changess = ["s1", "s2", "s1"]
    >>> changes_pos = [0, 3, 6, 10]
    >>> con_mat_col = {(0, 0): 0, (0, 1): 1, (0, 2): 2, (0, 3): 3, (0, 6): 6, (0, 7): 7, (0, 8): 8, (0, 9): 9, (1, 3): 1, (1, 4): 2, (1, 5): 3, (1, 6): 4}
    >>> strand_mat_col = {(0, 0): '+', (0, 1): '+', (0, 2): '+', (0, 3): '+', (0, 6): '+', (0, 7): '+', (0, 8): '+', (0, 9): '+', (1, 3): '-', (1, 4): '-', (1, 5): '-', (1, 6): '-'}
    >>> node_conf = {('s1', 0): 0.9, ('s1', 1): 0.5, ('s1', 2): 0.9, ('s1', 0): 0.1, ('s2', 1): 0.5, ('s2', 2): 0.1}
    >>> sub_col_ind = {'s1': 0, 's2': 1}
    >>> FillPathGraph(3, non_cols, changess, changes_pos, con_mat_col, strand_mat_col, node_conf, sub_col_ind)
    [0, 1, 1, 0, 0, 1, 0, 0, 0]
    """

    time1: float = time.time()

    path_graph: List[int] = []

    for i in range(nodes * nodes):
        path_graph.append(0)

    # filling beginning path graph with straight line through the nodes
    for i in range(nodes - 1):
        path_graph[i * nodes + i + 1] = 1

    for sink_node_index in range(2, nodes):  # don't need to add edges between first 2 nodes because already there
        sink_subfam: str = str(changes[sink_node_index])
        sink_subfam_start: int = consensus_matrix_collapse[subfams_collapse_index[sink_subfam], columns[changes_position[sink_node_index]]]
        sink_strand: str = strand_matrix_collapse[subfams_collapse_index[sink_subfam], columns[changes_position[sink_node_index]]]

        if sink_subfam != "skip":  # don't want to add alternative edges to skip nodes
            for source_node_index in range(sink_node_index - 1):
                source_subfam: str = changes[source_node_index]
                sourceConf: float = node_confidence[
                    sink_subfam, source_node_index]  # sink subfam confidence in source node
                sinkConf: float = node_confidence[
                    source_subfam, sink_node_index]  # source subfam confidence in sink node

                if (subfams_collapse_index[source_subfam], columns[changes_position[source_node_index + 1] - 1]) in consensus_matrix_collapse:

                    source_subfam_stop = consensus_matrix_collapse[
                        subfams_collapse_index[source_subfam], columns[changes_position[source_node_index + 1] - 1]]
                    source_strand = strand_matrix_collapse[
                        subfams_collapse_index[source_subfam], columns[changes_position[source_node_index + 1] - 1]]

                    # adds in edge if the subfam of the sink is at the source node and if it's
                    # confidence >= 30%, and if the source is before the sink in the consensus sequence

                    # FIXME - not sure what this confidence threshold should be
                    if sourceConf >= 0.2 or sinkConf >= 0.2:
                        # print(source_subfam, sink_subfam)
                        if sink_strand == '+' and sink_strand == source_strand:
                            # FIXME- not sure what this overlap should be .. just allowed 50 for now
                            if source_subfam_stop <= sink_subfam_start + 50:
                                path_graph[source_node_index * nodes + sink_node_index] = 1
                        elif sink_strand == '-' and sink_strand == source_strand:
                            if source_subfam_stop + 50 >= sink_subfam_start:
                                path_graph[source_node_index * nodes + sink_node_index] = 1
    print("FillPathGraph", time.time()- time1)

    return path_graph


def ExtractNodes(num_col: int, nodes: int, columns: List[int], changes_position: List[int],
                 path_graph: List[int]) -> int:
    """
    finds nodes that only have one (or less) incoming and one (or less) outgoing edge and adds
    them to RemoveStarts and RemoveStops so they can be extracted from the alignment - all columns
    in removed nodes are removed from NonEmptyColumns (columns) so during next round of DP calculations
    those nodes are ignored and surrounding nodes can be stitched if necessary

    input:
    num_col: number of columns in matrices
    nodes: number of nodes in graph
    columns: non empty columns in matrices
    changes_position: node boundaries
    path_graph: 2D array that represents edges in the graph

    output:
    updates columns
    updated_num_col: updated number of columns in matrices after nodes have been removed

    >>> non_cols = [0, 1, 2, 3, 4, 5, 6, 7, 8, 9]
    >>> change_pos = [0, 3, 7]
    >>> path_graph = [0, 1, 1, 0, 0, 1, 0, 0, 0]
    >>> ExtractNodes(10, 3, non_cols, change_pos, path_graph)
    10
    >>> non_cols
    [0, 1, 2, 7, 8, 9]
    """

    time1: float = time.time()

    remove_starts: List[int] = []
    remove_stops: List[int] = []

    updated_num_col: int = num_col

    # boolean for which nodes will be removed
    remove_nodes: List[bool] = [False for _ in range(nodes)]

    # extracting nodes that only have one incoming and one outgoing edge
    num_edges_in: List[int] = [0 for _ in range(nodes)]
    num_edges_out: List[int] = [0 for _ in range(nodes)]

    for row in range(nodes):
        for col in range(nodes):
            num_edges_in[col] += path_graph[row * nodes + col]
            num_edges_out[row] += path_graph[row * nodes + col]

    for node in range(nodes - 1):
        if num_edges_in[node] <= 1 and num_edges_out[node] <= 1:
            remove_starts.append(changes_position[node])
            remove_stops.append(changes_position[node + 1])
            remove_nodes[node] = True

    # deals with last node, so when NumNodes-1 the last remove stop is the end of the matrix
    if num_edges_in[nodes - 1] <= 1 and num_edges_out[nodes - 1] <= 1:
        remove_starts.append(changes_position[nodes - 1])
        remove_stops.append(len(columns) - 1)
        remove_nodes[nodes - 1] = True

    # when removing from the end, have to update num_cols because don't want the end of the matrix anymore
    col_index: int = nodes - 1
    while remove_nodes[col_index]:
        updated_num_col = columns[changes_position[col_index] - 1]
        col_index -= 1

    # removing inserted elements from NonEmptyColumns so they can be ignored
    total: int = 0
    for i in range(len(remove_stops)):
        del columns[remove_starts[i] - total:remove_stops[i] - total]
        # 	helps with offset, when first part is spliced out need an offset to know where to splice out for second part
        total += (remove_stops[i] - remove_starts[i])

    print("ExtractNodes", time.time() - time1)
    return updated_num_col


def PrintResults(changes_orig: List[str], changespos_orig: List[int], columns_orig: List[int], ids: List[int]) -> None:
    """
    prints the final results
    prints start and stop in terms of position in matrix
    """
    stdout.write("start\tstop\tID\tname\n")
    stdout.write("----------------------------------------\n")
    for i in range(len(changes_orig)):
        if str(changes_orig[i]) != 'skip':
            stdout.write(str(columns_orig[changespos_orig[i]]))
            stdout.write("\t")
            stdout.write(str(columns_orig[changespos_orig[i + 1] - 1]))
            stdout.write("\t")
            stdout.write(str(ids[columns_orig[changespos_orig[i]]]))
            stdout.write("\t")
            stdout.write(str(changes_orig[i]))
            stdout.write("\n")


def PrintResultsSequence(edgestart: int, changes_orig: List[str], changespos_orig: List[int], columns_orig: List[int],
                         ids: List[int]) -> None:
    """
    prints final results
    prints start and stop in terms of input chrom sequence
    """
    stdout.write("start\tstop\tID\tname\n")
    stdout.write("----------------------------------------\n")
    for i in range(len(changes_orig)):
        if str(changes_orig[i]) != 'skip':
            stdout.write(str(columns_orig[changespos_orig[i]] + edgestart))
            stdout.write("\t")
            stdout.write(str(columns_orig[changespos_orig[i + 1] - 1] + edgestart))
            stdout.write("\t")
            stdout.write(str(ids[columns_orig[changespos_orig[i]]]))
            stdout.write("\t")
            stdout.write(str(changes_orig[i]))
            stdout.write("\n")


def PrintResultsChrom(edgestart: int, chrom_start: int, changes_orig: List[str], changespos_orig: List[int], columns_orig: List[int],
                         ids: List[int]) -> None:
    """
    prints final results
    prints start and stop in terms of chromosome
    """
    stdout.write("start\tstop\tID\tname\n")
    stdout.write("----------------------------------------\n")
    for i in range(len(changes_orig)):
        if str(changes_orig[i]) != 'skip':
            stdout.write(str(columns_orig[changespos_orig[i]] + edgestart + chrom_start))
            stdout.write("\t")
            stdout.write(str(columns_orig[changespos_orig[i + 1] - 1] + edgestart + chrom_start))
            stdout.write("\t")
            stdout.write(str(ids[columns_orig[changespos_orig[i]]]))
            stdout.write("\t")
            stdout.write(str(changes_orig[i]))
            stdout.write("\n")


def PrintResultsViz(start_all: int, outfile: str, outfile_json: str, chrom: str, chrom_start: int, changes_orig: List[str], changes_position_orig: List[int],
                    columns_orig: List[int], consensus_lengths: Dict[str, int],
                    strand_matrix_collapse: Dict[Tuple[int, int], str],
                    consensus_matrix_collapse: Dict[Tuple[int, int], int], subfams_collapse_index: Dict[str, int]):
    """
    prints the results in the proper format to input into the SODA visualization tool

    format described here:
    https://genome.ucsc.edu/cgi-bin/hgTables?db=hg38&hgta_group=rep&hgta_track=joinedRmsk&hgta_table=rmskJoinedCurrent&hgta_doSchema=describe+table+schema

    """

    id: int = 0

    length = len(changes_position_orig) - 1

    used: List[int] = [1] * length  # wont print out the results of the same thing twice

    json_dict_id: Dict[str, Dict[str, Dict[str, float]]] = {}

    with open(outfile, 'w') as out:

        i = 0
        while i < length:

            sub_id: int = 0
            json_dict_subid: Dict[str, Dict[str, float]] = {}

            if changes_orig[i] != 'skip' and used[i]:
                subfam: str = changes_orig[i]
                strand: str = strand_matrix_collapse[subfams_collapse_index[subfam], columns_orig[changes_position_orig[i]]]
                consensus_start: int
                consensus_stop: int

                left_flank: int
                right_flank: int

                if strand == "-":
                    left_flank = consensus_lengths[subfam] - consensus_matrix_collapse[
                        subfams_collapse_index[subfam], columns_orig[changes_position_orig[i]]]
                    right_flank = consensus_matrix_collapse[
                                      subfams_collapse_index[subfam], columns_orig[changes_position_orig[i + 1] - 1]] - 1
                else:
                    left_flank = consensus_matrix_collapse[subfams_collapse_index[subfam], columns_orig[changes_position_orig[i]]] - 1
                    right_flank = consensus_lengths[subfam] - consensus_matrix_collapse[
                        subfams_collapse_index[subfam], columns_orig[changes_position_orig[i + 1] - 1]]

                align_start: int = chrom_start + (columns_orig[changes_position_orig[i]] + start_all)
                feature_start: int = align_start - left_flank
                align_stop: int = chrom_start + (columns_orig[changes_position_orig[i + 1] - 1] + start_all)
                feature_stop: int = align_stop + right_flank

                block_start_matrix: int = columns_orig[changes_position_orig[i]]

                block_count: int = 3
                id += 1

                json_dict_subfam_i: Dict[str, float] = {}

                for subfam_i in range(1, len(Subfams)):
                    subfamm = Subfams[subfam_i]
                    if NodeConfidenceOrig[subfamm, i] > 0.001:
                        json_dict_subfam_i[subfamm] = NodeConfidenceOrig[subfamm, i]

                json_dict_subid[str(id) + "-" + str(sub_id)] = sorted(json_dict_subfam_i.items(), key=lambda x: x[1], reverse=True)
                sub_id += 1

                block_start: List[str] = []
                block_size: List[str] = []

                block_start.append("-1")
                block_start.append(str(left_flank + 1))
                block_start.append("-1")

                block_size.append(str(left_flank))
                block_size.append(
                    str(columns_orig[changes_position_orig[i + 1] - 1] - columns_orig[changes_position_orig[i]] + 1))
                block_size.append(str(right_flank))

                j: int = i + 1
                while j < length:
                    if changes_orig[j] != 'skip':

                        if IDs[columns_orig[changes_position_orig[i]]] == IDs[columns_orig[changes_position_orig[j]]]:

                            if strand == "-":
                                right_flank = consensus_matrix_collapse[
                                                  subfams_collapse_index[changes_orig[j]], columns_orig[changes_position_orig[j + 1] - 1]] - 1
                            else:
                                right_flank = consensus_lengths[subfam] - consensus_matrix_collapse[
                                    subfams_collapse_index[subfam], columns_orig[changes_position_orig[j + 1] - 1]]

                            del block_size[-1]
                            block_size.append("0")

                            align_stop: int = chrom_start + (columns_orig[changes_position_orig[j + 1] - 1] + start_all)
                            feature_stop: int = align_stop + right_flank

                            block_start.append(
                                str(columns_orig[changes_position_orig[j]] + 1 - block_start_matrix + left_flank))
                            block_start.append("-1")

                            block_size.append(str(columns_orig[changes_position_orig[j + 1] - 1] - columns_orig[
                                changes_position_orig[j]] + 1))

                            block_size.append(str(right_flank))

                            block_count += 2

                            used[j] = 0

                            json_dict_subfam_j: Dict[str, float] = {}

                            for subfam_i in range(1, len(Subfams)):
                                subfamm = Subfams[subfam_i]
                                if NodeConfidenceOrig[subfamm, j] > 0.001:
                                    json_dict_subfam_j[subfamm] = NodeConfidenceOrig[subfamm, j]

                            json_dict_subid[str(id) + "-" + str(sub_id)] = sorted(json_dict_subfam_j.items(), key=lambda x: x[1], reverse=True)
                            sub_id += 1

                    j += 1

                json_dict_id[str(id)] = json_dict_subid

                out.write("000 " + chrom + " " + str(feature_start) + " " + str(
                    feature_stop) + " " + subfam + " 0 " + strand + " " + str(align_start) + " " + str(
                    align_stop) + " 0 " + str(block_count) + " " + (",".join(block_size)) + " " + (
                              ",".join(block_start)) + " " + str(id))
                out.write("\n")

            used[i] = 0
            i += 1

    #prints json file with confidence values for each annotation
    with open(outfile_json, 'w') as out_json:
        out_json.write(json.dumps(json_dict_id))



# -----------------------------------------------------------------------------------#
#			   MAIN														   			#
# -----------------------------------------------------------------------------------#


if __name__ == "__main__":

    # time1: float = time.time()
    time2: float = time.time()

    GapInit: int = -25
    GapExt: int = -5
    Lamb: float = 0.0
    EslPath = ""
    ChunkSize: int = 31
    SameProbLog: float = 0.0
    ChangeProb: float = 10 ** -45
    ChangeProbLog: float = 0.0  # Reassigned later
    ChangeProbSkip: float = 0.0  # Reassigned later
    SameProbSkip: float = 0.0
    SkipAlignScore: float = 0.0

    StartAll: int = 0  # Reassigned later
    StopAll: int = 0  # Reassigned later
    ID: int = 1111

    infile_prior_counts: str = ""
    outfile_viz: str = ""
    outfile_heatmap: str = ""

    help: bool = False  # Reassigned later
    prin: bool = False  # Reassigned later
    printMatrixPos: bool = False  # Reassigned later
    printSeqPos: bool = False

    helpMessage: str = f"""
    usage: {argv[0]} alignFile subMatrixFile\n
    ARGUMENTS
        --GapInit[-25]
        --getExt[-5]
        --lambda [will calc from matrix if not included]
        --eslPath [specify path to easel]
        --segmentsize[30]
        --changeprob[1e-45]
        --priorCounts PriorCountsFile
    
    OPTIONS
        --help - display help message
        --printmatrices - output all matrices
        --matrixpos - prints subfam changes in matrix position instead of genomic position
        --sequencepos - prints subfam changes in sequence position instead of genomic position
        --viz outfile - prints output format for SODA visualization
        --heatmap outfile - prints probability file for input into heatmap
    """

    raw_opts, args = getopt(argv[1:], "", [
        "GapInit=",
        "GapExt=",
        "skipScore=",
        "lambda=",
        "eslPath=",
        "segmentsize=",
        "changeprob=",
        "priorCounts=",
        "viz=",
        "heatmap=",

        "help",
        "matrixpos",
        "seqpos",
    ])
    opts = dict(raw_opts)

    GapInit = int(opts["--GapInit"]) if "--GapInit" in opts else GapInit
    GapExt = int(opts["--GapExt"]) if "--GapExt" in opts else GapExt
    Lamb = float(opts["--lambda"]) if "--lambda" in opts else Lamb
    EslPath = str(opts["--eslPath"]) if "--eslPath" in opts else EslPath
    ChunkSize = int(opts["--segmentsize"]) if "--segmentsize" in opts else ChunkSize
    ChangeProb = float(opts["--changeprob"]) if "--changeprob" in opts else ChangeProb
    infile_prior_counts = str(opts["--priorCounts"]) if "--priorCounts" in opts else infile_prior_counts
    outfile_viz = str(opts["--viz"]) if "--viz" in opts else outfile_viz
    outfile_heatmap = str(opts["--heatmap"]) if "--heatmap" in opts else outfile_heatmap
    help = "--help" in opts
    printMatrixPos = "--matrixpos" in opts
    printSeqPos = "--seqpos" in opts

    outfile_conf = outfile_viz + ".json"

    if help:
        print(helpMessage)
        exit(0)

    # input is alignment file of hits region and substitution matrix
    infile: str = args[0]
    infile_matrix: str = args[1]

    # Other open was moved down to where we load the alignments file
    with open(infile_matrix) as _infile_matrix:
        in_matrix: List[str] = _infile_matrix.readlines()

    if infile_prior_counts:
        with open(infile_prior_counts) as _infile_prior_counts:
            in_counts: List[str] = _infile_prior_counts.readlines()

    #if lambda isn't included at command line, run esl_scorematrix to calculate it from scorematrix
    if not Lamb:
        esl_stream = os.popen(EslPath + 'esl_scorematrix --dna 25p41g_edited.matrix')
        esl_output = esl_stream.read()
        esl_output_list = re.split(r"\n+", esl_output)
        lambda_list = re.split(r"\s+", esl_output_list[1])
        Lamb = float(lambda_list[2])

    # reads in the score matrix from file and stores in dict that maps 'char1char2' to the score from the
    # input substitution matrix - ex: 'AA' = 8
    SubMatrix: Dict[str, int] = {}
    line = in_matrix[0]
    line = re.sub(r"^\s+", "", line)
    line = re.sub(r"\s+$", "", line)
    chars = re.split(r"\s+", line)

    count: int = 0
    for line in in_matrix[1:]:
        line = re.sub(r"^\s+", "", line)
        line = re.sub(r"\s+$", "", line)
        subScores = re.split(r"\s+", line)
        for i in range(len(subScores)):
            SubMatrix[chars[count] + chars[i]] = int(subScores[i])
        count += 1

    # maps subfam names to genomic prior_count/total_in_genome from input file
    # used during confidence calculations
    SubfamCounts: Dict[str, float] = {}
    PriorTotal: float = 0
    prob_skip = 0.4  # about 60% of genome is TE derived #FIXME - do we want to hard code this in?
    if infile_prior_counts:
        for line in in_counts[1:]:
            line = re.sub(r"\n", "", line)
            info = re.split(r"\s+", line)
            SubfamCounts[info[0]] = int(info[1])
            PriorTotal += float(info[1])

        for key in SubfamCounts:
            SubfamCounts[key] = (1 - prob_skip) * SubfamCounts[key] / PriorTotal

        SubfamCounts["skip"] = prob_skip

    Subfams: List[str] = []
    Chroms: List[str] = []
    Scores: List[int] = []
    Strands: List[str] = []
    Starts: List[int] = []
    Stops: List[int] = []
    ConsensusStarts: List[int] = []
    ConsensusStops: List[int] = []
    SubfamSeqs: List[str] = []
    ChromSeqs: List[str] = []
    Flanks: List[int] = []

    AlignMatrix: Dict[Tuple[int, int], float] = {}
    # AlignMatrix: List[float] = []
    ConfidenceMatrix: Dict[Tuple[int, int], float] = {}
    SupportMatrix: Dict[Tuple[int, int], float] = {}
    ProbMatrix: Dict[Tuple[int, int], float] = {}
    OriginMatrix: Dict[Tuple[int, int], int] = {}
    ConsensusMatrix: Dict[Tuple[int, int], int] = {}
    SameSubfamChangeMatrix: Dict[Tuple[int, int], int] = {}

    NonEmptyColumns: List[int] = []
    ActiveCells: Dict[int, List[int]] = {}

    Changes: List[str] = []
    ChangesPosition: List[int] = []

    IDs: List[int] = []
    ChangesOrig: List[str] = []
    ChangesPositionOrig: List[int] = []
    NonEmptyColumnsOrig: List[int] = []

    SupportMatrixCollapse: Dict[Tuple[int, int], int] = {}
    ActiveCellsCollapse: Dict[int, List[int]] = {}
    SubfamsCollapse: List[str] = []
    SubfamsCollapseIndex: Dict[str, int] = {}
    ConsensusMatrixCollapse: Dict[Tuple[int, int], int] = {}
    StrandMatrixCollapse: Dict[Tuple[int, int], str] = {}

    # for graph/node part
    NumNodes: int = 0
    NodeConfidence: Dict[Tuple[str, int], float] = {}
    NodeConfidenceOrig: Dict[Tuple[str, int], float] = {}
    PathGraph: List[int] = []
    total: int = 0
    loop: int = 1

    # opens alignment file and stores all Subfams, Scores, Starts, Stops, subfams seqs and Chrom seqs in arrays
    numseqs: int = 0
    with open(infile) as _infile:
        alignments = load_alignments(_infile)
        for alignment in alignments:
            numseqs += 1

            Subfams.append(alignment.subfamily)
            Chroms.append(alignment.chrom)
            Scores.append(alignment.score)
            Strands.append(alignment.strand)
            Starts.append(alignment.start)
            Stops.append(alignment.stop)
            ConsensusStarts.append(alignment.consensus_start)
            ConsensusStops.append(alignment.consensus_stop)
            SubfamSeqs.append(alignment.subfamily_sequence)
            ChromSeqs.append(alignment.sequence)
            Flanks.append(alignment.flank)

    # if there is only one subfam in the alignment file, no need to run anything because we know
    # that subfam is what's there
    # numseqs = 2 because of the skip state
    if numseqs == 2:
        if printMatrixPos:
            stdout.write("start\tstop\tID\tname\n")
            stdout.write("----------------------------------------\n")
            stdout.write(f"{0}\t{Stops[1] - Starts[1]}\t1111\t{Subfams[1]}\n")
        else:
            stdout.write("start\tstop\tID\tname\n")
            stdout.write("----------------------------------------\n")
            stdout.write(f"{Starts[1]}\t{Stops[1]}\t1111\t{Subfams[1]}\n")
        exit()

    match = re.search(r"(.+):(\d+)-.+", Chroms[1])
    Chrom: str = match.groups()[0]
    ChromStart: int = int(match.groups()[1])

    ChangeProbLog = log(ChangeProb / (numseqs - 1))
    ChangeProbSkip = (ChangeProbLog / 2)  # jumping in and then out of the skip state counts as 1 jump
    SameProbSkip = ChangeProbLog / 20  # 5% of the jump penalty, staying in skip state for 20nt "counts" as one jump

    # precomputes number of rows and cols in matrices
    rows: int = len(Subfams)
    cols: int = 0  # assign cols in FillAlignMatrix

    # precomputes consensus seq length for PrintResultsViz()
    ConsensusLengths: Dict[str, int] = {}
    if outfile_viz:
        for i in range(1, len(Flanks)):
            if Strands[i] == "+":
                ConsensusLengths[Subfams[i]] = ConsensusStops[i] + Flanks[i]
            else:
                ConsensusLengths[Subfams[i]] = ConsensusStarts[i] + Flanks[i]

    # print("initialization and reading in files", time.time()-time1)
    # time1 = time.time()

    (StartAll, Stopall) = PadSeqs(Starts, Stops, SubfamSeqs, ChromSeqs)

    (cols, AlignMatrix) = FillAlignMatrix(StartAll, ChunkSize, GapExt, GapInit, SkipAlignScore, SubfamSeqs,
                                          ChromSeqs, Starts, SubMatrix)

    ConsensusMatrix = FillConsensusPositionMatrix(cols, rows, StartAll, SubfamSeqs, ChromSeqs, Starts, Stops, ConsensusStarts, Strands)

    (NonEmptyColumns, ActiveCells) = FillColumns(cols, rows, AlignMatrix)

    ConfidenceMatrix = FillConfidenceMatrix(rows, Lamb, infile_prior_counts, NonEmptyColumns, SubfamCounts, Subfams, ActiveCells,
                                            AlignMatrix)

    SupportMatrix = FillSupportMatrix(rows, cols, ChunkSize, StartAll, NonEmptyColumns, Starts, Stops, ConfidenceMatrix)

    (rows, ConsensusMatrixCollapse, StrandMatrixCollapse, SupportMatrixCollapse, SubfamsCollapse,
     ActiveCellsCollapse, SubfamsCollapseIndex) = CollapseMatrices(rows, NonEmptyColumns, Subfams, Strands, ActiveCells, SupportMatrix, ConsensusMatrix)

    if outfile_heatmap:
        PrintMatrixSupport(cols, StartAll, ChromStart, outfile_heatmap, SupportMatrixCollapse, SubfamsCollapse)

    (ProbMatrix, OriginMatrix, SameSubfamChangeMatrix) = FillProbabilityMatrix(SameProbSkip, SameProbLog, ChangeProbLog,
                                                                               ChangeProbSkip,
                                                                               NonEmptyColumns, SubfamsCollapse,
                                                                               ActiveCellsCollapse,
                                                                               SupportMatrixCollapse,
                                                                               StrandMatrixCollapse,
                                                                               ConsensusMatrixCollapse)

    IDs = [0] * cols

    (ID, ChangesPosition, Changes) = GetPath(cols, ID, NonEmptyColumns, IDs, ChangesOrig, ChangesPositionOrig,
                                             NonEmptyColumnsOrig, SubfamsCollapse, ActiveCellsCollapse, ProbMatrix,
                                             OriginMatrix, SameSubfamChangeMatrix)

    # keep the original annotation for reporting results
    ChangesOrig = Changes.copy()
    ChangesPositionOrig = ChangesPosition.copy()
    NonEmptyColumnsOrig = NonEmptyColumns.copy()

    # Steps-
    # 1.create confidence for nodes
    #	will be in a matrix that is #subfams x #nodes
    # 2.create path graph
    # 3.find alternative paths through graph and add those edges
    # 4.extract all nodes (from dp matrix) that have a single incoming and a single outgoing edge
    # 5.annotate again with removed subfams
    #   --stop when all nodes have incoming and outgoing edges <= 1 or there are <= 2 nodes left

    # print("original DP calculations", time.time() - time1)
    # time1 = time.time()

    count: int = 0
    while (True):
        count += 1
        NumNodes = len(Changes)

        # breakout of loop if there are 2 or less nodes left
        if (NumNodes <= 2):
            break

        NodeConfidence.clear()

        # initializes and fills node confidence matrix
        NodeConfidence = FillNodeConfidence(NumNodes, StartAll, GapInit, GapExt, Lamb, infile_prior_counts, NonEmptyColumns, Starts, Stops, ChangesPosition, Subfams, SubfamSeqs, ChromSeqs, SubfamCounts, SubMatrix)

        if count == 1:
            NodeConfidenceOrig = NodeConfidence.copy()

        PathGraph.clear()
        PathGraph = FillPathGraph(NumNodes, NonEmptyColumns, Changes, ChangesPosition,
                                  ConsensusMatrixCollapse, StrandMatrixCollapse, NodeConfidence, SubfamsCollapseIndex)

        # test to see if there nodes in the graph that have more that one incoming or outgoing edge,
        # if so keep looping, if not break out of the loop
        # if they are all 0, break out of the loop
        test: bool = False
        j: int = 0
        while j < NumNodes:
            i: int = 0
            while i < j - 1:
                if (PathGraph[i * NumNodes + j] == 1):
                    test = True
                i += 1
            j += 1

        if not test:
            break

        cols = ExtractNodes(cols, NumNodes, NonEmptyColumns, ChangesPosition, PathGraph)

        # using prob matrix and origin matrix, just skip the cols I'm not interested in and annotate
        # without the removed node - ignores removed nodes because they are not in NonEmptyColumns
        # using old prob matrix and origin matrix
        (ProbMatrix, OriginMatrix, SameSubfamChangeMatrix) = FillProbabilityMatrix(SameProbSkip, SameProbLog,
                                                                                   ChangeProbLog, ChangeProbSkip,
                                                                                   NonEmptyColumns, SubfamsCollapse,
                                                                                   ActiveCellsCollapse,
                                                                                   SupportMatrixCollapse,
                                                                                   StrandMatrixCollapse,
                                                                                   ConsensusMatrixCollapse)

        Changes.clear()
        ChangesPosition.clear()

        (ID, ChangesPosition, Changes) = GetPath(cols, ID, NonEmptyColumns, IDs, ChangesOrig, ChangesPositionOrig,
                                                 NonEmptyColumnsOrig, SubfamsCollapse, ActiveCellsCollapse, ProbMatrix,
                                                 OriginMatrix, SameSubfamChangeMatrix)

    # print("graph stitching", time.time() - time1)
    # time1 = time.time()

    if printMatrixPos:
        PrintResults(ChangesOrig, ChangesPositionOrig, NonEmptyColumnsOrig, IDs)
    elif printSeqPos:
        PrintResultsSequence(StartAll, ChangesOrig, ChangesPositionOrig, NonEmptyColumnsOrig, IDs)
    else:
        PrintResultsChrom(StartAll, ChromStart, ChangesOrig, ChangesPositionOrig, NonEmptyColumnsOrig, IDs)

    if outfile_viz:
        PrintResultsViz(StartAll, outfile_viz, outfile_conf, Chrom, ChromStart, ChangesOrig, ChangesPositionOrig, NonEmptyColumnsOrig, ConsensusLengths,
                        StrandMatrixCollapse, ConsensusMatrixCollapse, SubfamsCollapseIndex)

    # print("print results", time.time() - time1)

    print("ALL", time.time()-time2)

