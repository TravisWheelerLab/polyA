from typing import Dict
from math import log


def calculate_score(
    gap_ext: int,
    gap_init: int,
    seq1: str,
    seq2: str,
    prev_char_seq1: str,
    prev_char_seq2: str,
    sub_matrix: Dict[str, int],
) -> float:
    """
    Calculate the score for an alignment between a subfamily and a target/chromsome sequence.
    Scores are calculated based on input sub_matrix, gap_ext, and gap_init.

    prev_char_seq1, prev_char_seq2 are the single nucleotides in the alignment before the chunk - if
    chunk starts with '-' these tell us to use gap_init or gap_ext as the penalty

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
    >>> calculate_score(-5, -25, "AT", "AT", "", "", sub_mat)
    2.0
    >>> calculate_score(-5, -25, "-T", "AT", "A", "A", sub_mat)
    -24.0
    >>> calculate_score(-5, -25, "-T", "AT", "-", "", sub_mat)
    -4.0
    """
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
    else:
        chunk_score += sub_matrix[seq1[0] + seq2[0]]

    for j in range(1, len(seq1)):
        seq1_char: str = seq1[j]
        seq2_char: str = seq2[j]

        if seq1_char == "-":
            if seq1[j - 1] == "-":
                chunk_score += gap_ext
            else:
                chunk_score += gap_init
        elif seq2_char == "-":
            if seq2[j - 1] == "-":
                chunk_score += gap_ext
            else:
                chunk_score += gap_init
        else:
            chunk_score += sub_matrix[seq1_char + seq2_char]

    return float(chunk_score)


def calculate_hmm_score(
    hmm_start: int,
    chrom: str,
    subfam: str,
    prev_char_seq1: str,
    prev_char_seq2: str,
    subfam_hmm,
) -> float:
    """
    Function description here
    hmm start - hmm pos of first char
    """
    chunk_score: float = 0
    return float(chunk_score)


def calculate_insertion_score(
    hmm_pos: int,
    chrom: str,
    subfam: str,
    subfam_hmm,
) -> float:
    """
    Function description here

    """
    i: int = 0
    insertion_score: float = float(subfam_hmm[hmm_pos]["transition"]["m->i"])
    while i + 1 < len(subfam) and subfam[i + 1] == "-":
        insertion_score += float(subfam_hmm[hmm_pos]["transition"]["i->i"])
        i += 1
    i += 1
    insertion_score += float(subfam_hmm[hmm_pos]["transition"]["i->m"])
    return insertion_score / float(i)
