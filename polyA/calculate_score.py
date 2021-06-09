from typing import Dict, Optional, Tuple, Set, List
import math
from .constants import CROSS_MATCH_ADJUSTMENT


def calculate_score(
    gap_ext: float,
    gap_init: float,
    seq1: str,
    seq2: str,
    prev_char_seq1: str,
    prev_char_seq2: str,
    sub_matrix: Dict[str, int],
    char_complexity_adjustment: Dict[str, float],
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
    >>> char_adjustments = {'T': 0.0, 'C': 0.0, 'A': 0.0, 'G': 0.0, '-': 0.0}
    >>> calculate_score(-5, -25, "AT", "AT", "", "", sub_mat, char_adjustments)
    2.0
    >>> calculate_score(-5, -25, "-T", "AT", "A", "A", sub_mat, char_adjustments)
    -24.0
    >>> calculate_score(-5, -25, "-T", "AT", "-", "", sub_mat, char_adjustments)
    -4.0
    >>> char_adjustments = {'T': -0.4, 'C': 0.0, 'A': -0.2, 'G': 0.0, '-': 0.0}
    >>> calculate_score(-5, -25, "-T", "AT", "A", "A", sub_mat, char_adjustments)
    -24.4
    >>> calculate_score(-5, -25, "-T", "AT", "-", "", sub_mat, char_adjustments)
    -4.4
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
            chunk_score += char_complexity_adjustment[seq2_char]
    return chunk_score


def calc_target_char_counts(
    query_seq: str,
    target_seq: str,
    target_chars: List[str],
) -> Tuple[Dict[str, int], Set[str]]:
    """
    Finds the char counts in the target sequence which will be used in the calculation
    to determine a character's contribution to the complexity adjusted score.

    input:
    query_seq: query sequence from alignment
    target_seq: target sequence from alignment
    target_chars: chars in the target sequence that will contribute to
    the complexity adjusted score

    output:
    target_char_counts: dictionary of the target chars to their count
    in the target sequence that will contribute to the adjusted score
    other_chars: other chars found in the target sequence that do not contribute
    to the complexity adjustment (gaps, padding)
    """
    target_char_counts: Dict[str, int] = {char: 0 for char in target_chars}
    other_chars = set()
    for q_base, t_base in zip(query_seq, target_seq):
        if t_base not in target_chars:
            other_chars.add(t_base)
            continue
        if q_base == "-":
            continue
        target_char_counts[t_base] += 1
    return target_char_counts, other_chars


def calculate_complexity_adjusted_score(
    char_background_freqs: Optional[Dict[str, float]],
    query_seq: str,
    target_seq: str,
    lamb: float,
) -> Dict[str, float]:
    r"""
    Calculates the per position contribution to the complexity adjusted score
    of the alignment for each valid char in the target sequence.
    This calculation is from the cross_match complexity adjustment functions:

    .. math::
        t_{i} = \text{count of $i^{th}$ character}

    .. math::
        p_{i} = \text{background frequency of $i^{th}$ character}

    .. math::
        t_{factor} = \sum_{i=0}^{N} t_{i}ln(t_{i}) - \left[ \left( \sum_{i=0}^{N} t_{i} \right) ln(\sum_{i=0}^{N} t_{i}) \right]

    .. math::
        t_{sum} = \left[ \sum_{i=0}^{N} t_{i} ln(p_{i}) \right] - t_{factor}

    .. math::
        s_{adj} = s_{raw} + (\frac{t_{sum}}{\lambda}) + 0.999

    input:
    char_background_freqs: background frequencies of the chars from the scoring matrix
    query_seq: query sequence from alignment
    target_seq: target sequence from alignment
    lamb: lambda value used from the scoring matrix

    output:
    char_complexity_adjustments: a char's per position contribution to the
    complexity adjusted score

    >>> query = "TCAGACTGTTCA-----ACTCACCTGGCAGCCACTTCCAGA"
    >>> target = "TCAGACTGTTCATGAGTGCTCACCTGGTAGAGG-----AAA"
    >>> lamb = 0.1227
    >>> calculate_complexity_adjusted_score(None, query, target, lamb)
    {'T': 0.0, 'C': 0.0, 'A': 0.0, 'G': 0.0, '-': 0.0}

    >>> char_freqs = {'A': 0.295, 'G': 0.205, 'C': 0.205, 'T': 0.295}
    >>> calculate_complexity_adjusted_score(char_freqs, query, target, lamb)
    {'A': -0.08342236356495786, 'G': -0.11790190327726593, 'C': -0.11790190327726593, 'T': -0.08342236356495786, '-': 0.0}
    """
    char_complexity_adjustments: Dict[str, float] = {}
    if char_background_freqs is None:
        # set the complexity adjustment to zero for every char
        for char in target_seq:
            char_complexity_adjustments[char] = 0.0
        return char_complexity_adjustments

    t_factor = 0.0
    t_sum = 0.0
    t_counts: int = 0
    target_chars = list(char_background_freqs.keys())
    target_char_counts, other_chars = calc_target_char_counts(
        query_seq, target_seq, target_chars
    )
    total_chars = sum(target_char_counts.values())

    if total_chars > 0:
        char_complexity_adjustments = {char: 0 for char in target_chars}
        for char, freq in char_background_freqs.items():
            count = target_char_counts[char]
            if count > 0 and 0 < freq < 1:
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
                char_complexity_adjustments[c] += (
                    CROSS_MATCH_ADJUSTMENT / total_chars
                )

    # to avoid any key errors
    for char in other_chars:
        char_complexity_adjustments[char] = 0.0

    return char_complexity_adjustments
