from typing import Dict, Optional
import math


def calculate_score(
    gap_ext: float,
    gap_init: float,
    seq1: str,
    seq2: str,
    prev_char_seq1: str,
    prev_char_seq2: str,
    sub_matrix: Dict[str, int],
    char_complexity_adjustment: Dict[str, int] = None,
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
        if char_complexity_adjustment:
            chunk_score += char_complexity_adjustment[seq2[0]]

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
            if char_complexity_adjustment:
                chunk_score += char_complexity_adjustment[seq2_char]
    return chunk_score


def calc_target_char_counts(
    query_seq: str,
    target_seq: str,
):
    query_char_counts = {"A": 0, "C": 0, "G": 0, "T": 0}
    other_chars = set()
    total_chars = 0
    for q_base, t_base in zip(query_seq, target_seq):
        if t_base not in ["A", "G", "T", "C"]:
            other_chars.add(t_base)
            continue
        if q_base == "-":
            continue
        query_char_counts[t_base] += 1
        total_chars += 1
    return query_char_counts, total_chars, other_chars


def calculate_complexity_adjusted_score(
    char_background_freqs: Optional[Dict[str, float]],
    query_seq: str,
    target_seq: str,
    lamb: float,
):
    char_complexity_adjustments: Dict[str, float] = {}
    if char_background_freqs is None:
        # set the complexity adjustment to zero for every char
        for char in target_seq:
            char_complexity_adjustments[char] = 0
            if char == ".":
                break
        return char_complexity_adjustments

    t_factor: float = 0
    t_sum: float = 0
    t_counts: int = 0
    target_char_counts, total_chars, other_chars = calc_target_char_counts(
        query_seq, target_seq
    )

    char_complexity_adjustments = {"A": 0, "C": 0, "G": 0, "T": 0}
    for char, freq in char_background_freqs.items():
        count = target_char_counts[char]
        if count > 0 and freq > 0 and math.log(freq) != 0:
            t_factor += count * math.log(count)
            t_sum += count * math.log(freq)
            t_counts += count
            char_complexity_adjustments[char] += count * math.log(freq)

    # char percent contribution to t_sum
    for c, v in char_complexity_adjustments.items():
        char_complexity_adjustments[c] /= t_sum

    if t_counts != 0:
        t_factor -= t_counts * math.log(t_counts)
    t_sum -= t_factor

    # char value contribution to t_sum
    for c, v in char_complexity_adjustments.items():
        char_complexity_adjustments[c] *= t_sum

    # per position char value contribution to score adjustment
    for c, v in char_complexity_adjustments.items():
        if v != 0:
            char_complexity_adjustments[c] /= target_char_counts[c]
            char_complexity_adjustments[c] /= lamb
            char_complexity_adjustments[c] += 0.999 / total_chars

    # to avoid any key errors
    for char in other_chars:
        char_complexity_adjustments[char] = 0

    return char_complexity_adjustments
