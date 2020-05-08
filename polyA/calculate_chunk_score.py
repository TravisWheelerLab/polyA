import logging

from substitution_matrix import SubstitutionMatrix

_logger = logging.getLogger(__name__)


def calculate_chunk_score(
    subfamily_chunk: str,
    sequence_chunk: str,
    gap_extend_score: int,
    gap_start_score: int,  #TODO: GEORGE - change name to gap_init_score
    substitution_matrix: SubstitutionMatrix,
    prev_subfamily_chunk: str = "", #TODO: GEORGE - prev_subfamily_char
    prev_sequence_chunk: str = "",	#TODO: GEORGE - prev_sequence_char
) -> int:
    """
    Calculate the score for a particular alignment between a subfamily
    and a sequence. Scores are calculated based on input SubstitutionMatrix, gap_extend_score,
    and gap_init_score.
    
    prev_subfamily_char, prev_sequence_char are the single nucleotides in the alignment 
    before the chunk - if chunk starts with '-' these tell us to use gap_init_score or 
    gap_extend_score as the penalty

    >>> sub_mat = {("a", "a"): 1, ("a", "b"): 2, ("b", "a"): 4, ("b", "b"): 8}
    >>> calculate_chunk_score("ab", "ab", 0, 0, sub_mat)
    9
    >>> calculate_chunk_score("ba", "ab", 0, 0, sub_mat)
    6
    >>> calculate_chunk_score("-b", "ab", 0, 10, sub_mat, prev_subfamily_chunk = "a")
    18
    >>> calculate_chunk_score("-b", "ab", 10, 0, sub_mat, prev_subfamily_chunk = "-")
    18
    """
    chunk_score: int = 0

    # If the first character in the subfamily chunk is empty,
    # then we need to look back at the previous subfamily chunk
    # to determine whether this is the start of a new gap or an
    # extension of an old gap.
    if subfamily_chunk[0] == "-":
        if prev_subfamily_chunk[-1] == "-":
            chunk_score += gap_extend_score
        else:
            chunk_score += gap_start_score
    elif sequence_chunk[0] == "-":
        if prev_sequence_chunk[-1] == "-":
            chunk_score += gap_extend_score
        else:
            chunk_score += gap_start_score
    elif subfamily_chunk[0] == "." or sequence_chunk[0] == ".":
        # One or the other of the chunks begins with padding, so
        # we can't make any determination about the score.
        chunk_score = chunk_score
    else:
        chunk_score += substitution_matrix[
            subfamily_chunk[0], sequence_chunk[0]
        ]

    for j in range(1, len(subfamily_chunk)):
        if subfamily_chunk[j] == "-":
            if subfamily_chunk[j - 1] == "-":
                chunk_score += gap_extend_score
            else:
                chunk_score += gap_start_score
        elif sequence_chunk[j] == "-":
            if sequence_chunk[j - 1] == "-":
                chunk_score += gap_extend_score
            else:
                chunk_score += gap_start_score
        elif subfamily_chunk[j] == "." or sequence_chunk[j] == ".":
            chunk_score = chunk_score
        else:
            # TODO: Decide how to handle this situation and justify it
            nuc_pair = subfamily_chunk[j], sequence_chunk[j]
            if nuc_pair in substitution_matrix:
                chunk_score += substitution_matrix[nuc_pair]
            else:
                _logger.warning(
                    f"nucleotide pair {nuc_pair} missing from substitution matrix"
                )

    return chunk_score


if __name__ == "__main__":
    import doctest

    doctest.testmod()
