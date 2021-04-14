from typing import Dict
import math


def calculate_score(
    gap_ext: float,
    gap_init: float,
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
    chunk_score: float = 0.0

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

    return chunk_score


def calc_query_char_counts(query_seq, target_seq):
    query_char_count = {"A": 0, "C": 0, "G": 0, "T": 0}
    total_chars = 0
    for q_base, t_base in zip(target_seq, query_seq):
        if q_base not in ["A", "G", "T", "C"]:
            continue
        if t_base == "-":
            continue
        query_char_count[q_base] += 1
        total_chars += 1
    return query_char_count, total_chars


def calculate_complexity_adjusted_score(query_seq, target_seq, lmbda, score):
    t_factor = 0
    t_sum = 0
    t_counts = 0
    # background_freq from CM file
    background_freq = {"A": 0.295, "C": 0.205, "G": 0.205, "T": 0.295}
    query_char_counts, total_chars = calc_query_char_counts(
        query_seq, target_seq
    )
    position_char_adjustment = {"A": 0, "C": 0, "G": 0, "T": 0}
    for char, freq in background_freq.items():
        count = query_char_counts[char]
        if count > 0 and freq > 0 and math.log(freq) != 0:
            t_factor += count * math.log(count)
            t_sum += count * math.log(freq)
            t_counts += count
            position_char_adjustment[char] += count * math.log(freq)

    # char percent contribution to t_sum
    for c, v in position_char_adjustment.items():
        position_char_adjustment[c] /= t_sum

    if t_counts != 0:
        t_factor -= t_counts * math.log(t_counts)
    t_sum -= t_factor

    # char value contribution to t_sum
    for c, v in position_char_adjustment.items():
        position_char_adjustment[c] *= t_sum

    # per position char value contribution to score adjustment
    for c, v in position_char_adjustment.items():
        position_char_adjustment[c] /= query_char_counts[c]
        position_char_adjustment[c] /= lmbda
        position_char_adjustment[c] += 0.999 / total_chars

    # check that summ matches t_sum / lmbda + 0.999
    # FIXME: slightly off for shorter lengths - rounding errors?
    summ = 0
    for q, t in zip(query_seq, target_seq):
        if q not in ["A", "G", "T", "C"]:
            continue
        if t == "-":
            continue
        summ += position_char_adjustment[q]

    new_score = int(score + t_sum / lmbda + 0.999)
    new_score = max(0, new_score)
    # FIXME: do I need to check if 0 is returned?
    return position_char_adjustment
