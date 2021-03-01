from typing import List, Tuple


def edges(starts: List[int], stops: List[int]) -> Tuple[int, int]:
    """
    Find and return the min start and max stop positions for the entire region
    included in the alignments.

    Inputs:

    starts - start positions on the target sequence from the input alignment
    stops - stop positions on the target sequence from the input alignment

    Outputs:

    minimum and maximum start and stop positions on the chromosome/target
    sequences for whole alignment

    >>> b, e = edges(starts=[0, 1, 4, 7], stops=[0, 3, 10, 9])
    >>> b
    1
    >>> e
    10
    """
    min_start: int = min(starts[1:])
    max_stop: int = max(stops[1:])

    return min_start, max_stop
