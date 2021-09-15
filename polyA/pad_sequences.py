from typing import List


def pad_sequences(
    chunk_size: int,
    subfam_seqs: List[str],
    chrom_seqs: List[str],
    padding_char: str = ".",
) -> None:
    """
    Right-pad sequences with `(chunk_size - 1) / 2` copies of ".".

    The first sequence is ignored since that is the skip state.

    The `padding_char` parameter exists because the doctest doesn't "see" the
    periods, so the test will pass no matter how many are added to each
    sequence.

    Inputs:

    chunk_size: size (in nucleotides) of each segment
    subfam_seqs: actual subfamily / consensus sequences from alignment
    chrom_seqs: actual target / chromosome sequences from alignment
    padding_char: for testing, specifies the character to use for padding

    Side effects:

    Updates subfam_seqs and chrom_seqs with padded sequences

    >>> s_seq = ['', 'a', 'aaa']
    >>> c_seq = ['', 'a', 't-t']
    >>> pad_sequences(31, s_seq, c_seq, padding_char="-")
    >>> s_seq
    ['', 'a----------------', 'aaa----------------']
    >>> c_seq
    ['', 'a----------------', 't-t----------------']
    """
    padding = padding_char * ((chunk_size + 1) // 2)

    for i in range(1, len(subfam_seqs)):
        subfam_seqs[i] = subfam_seqs[i] + padding
        chrom_seqs[i] = chrom_seqs[i] + padding
