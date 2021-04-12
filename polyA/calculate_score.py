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


def query_char_frequency(char, query_seq):
    char_count = 0
    total_chars = 0
    for query_char in query_seq:
        if query_char == char:
            char_count += 1
        if query_char != "-":
            total_chars += 1
    # want the amount of that char / len(query_seq) excluding gaps
    return char_count / total_chars


def calculate_complexity_adjusted_score(query_seq, lmbda, score):
    t_factor = 0
    t_sum = 0
    t_counts = 0
    # From CM file
    background_freq = {"A": 0.295, "C": 0.205, "G": 0.205, "T": 0.295}
    for char, freq in background_freq.items():
        count = query_char_frequency(char, query_seq)
        if count > 0 and freq > 0 and math.log(freq) != 0:
            t_factor += count * math.log(count)
            t_sum += count * math.log(freq)
            t_counts += count

    if t_counts != 0:
        t_factor -= t_counts * math.log(t_counts)
    t_sum -= t_factor
    new_score = int(score + (t_sum / lmbda) + 0.999)
    new_score = max(0, new_score)
    return new_score
