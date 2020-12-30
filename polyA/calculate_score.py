from typing import Dict


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
    subfam: ful subfam seq
    insertion_score:
    insertion_index:
    subfam_hmm: dictionary of hmm family info

    output:
    HMM alignment score
    """
    chunk_score: float = 0
    hmm_pos = hmm_start
    for i in range(len(chrom_slice)):
        if chrom_slice[i] == ".":
            break
        elif chrom_slice[i] == "-":
            # deletion score
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
                insertion_score = calculate_insertion_score(
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


def calculate_insertion_score(
    hmm_pos: int,
    model: str,
    subfam_hmm,
) -> float:
    """
    Calculates the full insertion score in the model sequence

    input:
    hmm pos: hmm start pos of first nuc in model
    model: sequence - starts with the first gap
    subfam_hmm: dictionary of hmm family info

    output:
    per gap insertion score
    """
    i: int = 0  # num of gaps in model
    insertion_score: float = float(subfam_hmm[hmm_pos]["transition"]["m->i"])
    while model[i + 1] == "-":
        insertion_score += float(subfam_hmm[hmm_pos]["transition"]["i->i"])
        i += 1
    i += 1
    insertion_score += float(subfam_hmm[hmm_pos]["transition"]["i->m"])
    return insertion_score / float(i)
