from typing import Dict, Tuple


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
    chrom_slice: str,
    subfam_slice: str,
    subfam: str,
    insertion_score: float,
    insertion_index: int,
    subfam_hmm,
) -> float:
    """
    Calculate the score for an HMM alignment between a subfamily/model
    and a target/chromosome sequence.

    input:
    hmm_start: hmm start pos of first char in model
    chrom_slice:
    subfam_slice:
    subfam: full subfam seq starting at the same pos of subfam slice
    insertion_score: per gap score of full insertion
    insertion_index: index of last gap char relative to subfam slice
    subfam_hmm: dictionary of hmm family info

    output:
    HMM alignment score
    1112
    A--T

    1123
    --TT
    """
    chunk_score: float = 0
    hmm_pos = hmm_start
    for i in range(len(chrom_slice)):
        if chrom_slice[i] == ".":
            break
        elif chrom_slice[i] == "-":
            # deletion score
            # must have a gap or nuc in prev index
            if chrom_slice[i - 1] != "-":
                # match to deletion from prev hmm pos
                chunk_score += float(
                    subfam_hmm[hmm_pos - 1]["transition"]["m->d"]
                )
            else:
                # deletion to deletion from prev hmm pos
                chunk_score += float(
                    subfam_hmm[hmm_pos - 1]["transition"]["d->d"]
                )
            if chrom_slice[i + 1] != "-":
                # deletion to match from cur hmm pos
                chunk_score += float(subfam_hmm[hmm_pos]["transition"]["d->m"])
            hmm_pos += 1
        elif subfam_slice[i] == "-":
            # insertion score
            if insertion_index + 1 != i:
                # new insertion
                insertion_score = calculate_new_insertion_score(
                    hmm_pos, subfam[i::], subfam_hmm
                )
            insertion_index = i
            chunk_score += insertion_score
        else:
            # match
            chunk_score += float(
                subfam_hmm[hmm_pos]["emission"][chrom_slice[i]]
            )
            hmm_pos += 1
    return chunk_score


def calculate_new_insertion_score(
    hmm_pos: int,
    model_slice: str,
    subfam_hmm,
) -> float:
    """
    Calculates a new insertion score in the model sequence

    input:
    hmm pos: hmm start pos of first char in model
    model_slice: sequence - starts with the first gap, ex: GA[---AT....]
    subfam_hmm: dictionary of hmm family info

    output:
    per gap insertion score

    hmm pos: 223
    subfam:  --A
    m->i from 2, i->i from 2, i->i from 2, i->m from 2
    """
    i: int = 0  # num of gaps in model
    insertion_score: float = float(subfam_hmm[hmm_pos]["transition"]["m->i"])
    while model_slice[i + 1] == "-":  # will eventually end: --A
        insertion_score += float(subfam_hmm[hmm_pos]["transition"]["i->i"])
        i += 1
    i += 1
    insertion_score += float(subfam_hmm[hmm_pos]["transition"]["i->m"])
    return insertion_score / float(i)


def calculate_full_insertion_score(
    start_index: int,
    hmm_pos: int,
    model: str,
    subfam_hmm,
) -> Tuple[int, float]:
    """
    Calculates a full insertion score in the model
    starting at an index of a gap

    input:
    start_index: index in model of start of model slice
    hmm pos: hmm pos of first char in model_slice
    model: full subfam model sequence
    subfam_hmm: dictionary of hmm subfam info

    output:
    index of prev gap, per gap insertion score
    """
    # from current position (start_index) - look forwards and backwards
    # if gap at model_slice[first_index - 1] -> return first_index - 1
    # must be m->i and i->m
    insertion_score: float = float(subfam_hmm[hmm_pos]["transition"]["m->i"])
    gap_count: int = 1
    # search forward in model from start_index
    for i in range(start_index + 1, len(model)):
        if model[i] == "-":
            gap_count += 1
        else:
            break
    # search forward in model from start_index
    prev_gap: int = -2
    for i in range(start_index - 1, -1):
        if model[i] == "-":
            prev_gap = start_index - 1
            gap_count += 1
        else:
            break
    insertion_score += float(subfam_hmm[hmm_pos]["transition"]["i->i"]) * (
        gap_count - 1
    )
    insertion_score += float(subfam_hmm[hmm_pos]["transition"]["i->m"])
    return prev_gap, insertion_score / float(gap_count)
