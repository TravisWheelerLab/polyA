from typing import List, Tuple


def edges(starts: List[int], stops: List[int]) -> Tuple[int, int]:
    """
    Find and return the min start and max stop positions for the entire region
    included in the alignments.

    input:
    starts: start positions on the target sequence from the input alignment
    stops: stop positions on the target sequence from the input alignment

    output:
    minimum and maximum start and stop positions on the chromosome/target sequences for whole alignment

    >>> starts = [0, 1, 4, 7]
    >>> stops = [0, 3, 10, 9]
    >>> b, e = edges(starts, stops)
    >>> b
    1
    >>> e
    10
    """
    min_start: int = starts[1]
    max_stop: int = stops[1]

    for i in range(1, len(starts)):
        if starts[i] < min_start:
            min_start = starts[i]
        if stops[i] > max_stop:
            max_stop = stops[i]

    return min_start, max_stop
