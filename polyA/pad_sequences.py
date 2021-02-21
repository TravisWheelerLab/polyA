from typing import Dict, List, Tuple


def pad_sequences(
    chunk_size: int,
    subfam_seqs: List[str],
    chrom_seqs: List[str],
) -> None:
    """
    Pad out sequences with "." to allow regions where sequences do not all
    have the same start and stop positions.

    pad with an extra (chunk_size-1)/2 at the end

    input:
    chunk_size: size (in nucleotides) of each segment
    start: start positions on the target sequence from the input alignment
    stop: stop positions on the target sequence from the input alignment
    subfam_seqs: actual subfamily/consensus sequences from alignment
    chrom_seqs: actual target/chromosome sequences from alignment

    output:
    updates subfam_seqs and chrom_seqs with padded sequences
    minimum and maximum start and stop positions on the chromosome/target sequences for whole alignment

    >>> s_seq = ['', 'a', 'aaa']
    >>> c_seq = ['', 'a', 't-t']
    >>> pad_sequences(31, s_seq, c_seq)
    >>> s_seq
    ['', 'a.....................................', 'aaa.................................']
    >>> c_seq
    ['', 'a.....................................', 't-t.................................']
    """
    for i in range(1, len(subfam_seqs)):
        chrom_seqs[i] = f"{chrom_seqs[i]}" + ("." * chunk_size)
        subfam_seqs[i] = f"{subfam_seqs[i]}" + ("." * chunk_size)
