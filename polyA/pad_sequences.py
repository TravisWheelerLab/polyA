from typing import List, Tuple

from polyA.edges import edges


def pad_sequences(
    chunk_size: int,
    start: List[int],
    stop: List[int],
    subfam_seqs: List[str],
    chrom_seqs: List[str],
) -> Tuple[int, int]:
    """
    Pad out sequences with "." to allow regions where sequences do not all
    have the same start and stop positions.

    pad with an extra (chunk_size-1)/2 at the end

    input:
    start: start positions on the target sequence from the input alignment
    stop: stop positions on the target sequence from the input alignment
    subfam_seqs: actual subfamily/consensus sequences from alignment
    chrom_seqs: actual target/chromosome sequences from alignment

    output:
    updates subfam_seqs and chrom_seqs with padded sequences
    minimum and maximum start and stop positions on the chromosome/target sequences for whole alignment

    >>> starts = [0, 1, 3]
    >>> stops = [0, 1, 5]
    >>> s_seq = ['', 'a', 'aaa']
    >>> c_seq = ['', 'a', 't-t']
    >>> (b, e) = pad_sequences(31, starts, stops, s_seq, c_seq)
    >>> b
    1
    >>> e
    5
    >>> s_seq
    ['', 'a...................', '..aaa...............']
    >>> c_seq
    ['', 'a...................', '..t-t...............']
    """
    edge_start: int
    edge_stop: int

    half_chunk: int = int((chunk_size - 1) / 2)

    edge_start, edge_stop = edges(start, stop)

    for i in range(1, len(subfam_seqs)):
        left_pad: int = start[i] - edge_start
        right_pad: int = edge_stop - stop[i]

        chrom_seqs[i] = (
            ("." * left_pad)
            + f"{chrom_seqs[i]}"
            + ("." * (right_pad + half_chunk))
        )
        subfam_seqs[i] = (
            ("." * left_pad)
            + f"{subfam_seqs[i]}"
            + ("." * (right_pad + half_chunk))
        )

    return edge_start, edge_stop
