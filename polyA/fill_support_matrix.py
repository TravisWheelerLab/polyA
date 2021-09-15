from typing import Dict, List, Tuple

from .matrices import ConfidenceMatrix, SupportMatrix
from .performance import timeit


@timeit()
def fill_support_matrix(
    chunk_size: int,
    starts: List[int],
    stops: List[int],
    confidence_matrix: ConfidenceMatrix,
) -> SupportMatrix:
    """
    Fill the support_matrix using values in `confidence_matrix`. Average
    confidence values for surrounding `chunk_size` confidence values -
    normalized by dividing by number of segments.

    The score for subfam x at position i is the sum of all confidences for
    subfam x for all segments that overlap position i - divided by the number of
    segments.

    Inputs:

    chunk_size - the width of the window to use when averaging confidence
    scores. This is assumed to be an odd number >= 3.

    starts - a list of start columns for each row in the confidence matrix.
    These positions come should come directly from the alignments without any
    offsetting. The exception is the first value, which should contain the
    minimum of the remaining values, less one, in other words, `starts[0] ==
    min(starts[1:]) - 1`. This is the implicit start position for the skip
    state, which occupies the first row and is the only row that has data in the
    first column of the confidence matrix.

    stops - equivalent to the starts list, but contains the last column in each
    row that contains data. Again, the first value is special and should contain
    the maximum of the remaining values. In other words, `stops[0] ==
    max(stops[1:]) + 1`.

    confidence_matrix - confidence values computed with fill_confidence_matrix.
    The confidence matrix should have one contiguous run of values per row.

    Outputs:

    support_matrix - support values represent average confidence values over a
    window chunk_size wide at each position in the confidence matrix. This is
    the "support score".

    >>> conf_mat = {
    ...     (0, 0): 1.0, (0, 1): 0.1, (0, 2): 0.2, (0, 3): 0.4, (0, 4): 0.5, (0, 5): 1.0,
    ...     (1, 1): 0.9, (1, 2): 0.5, (1, 3): 0.5, (1, 4): 0.5,
    ...     (2, 2): 0.3, (2, 3): 0.1,
    ... }
    >>> supp_mat = fill_support_matrix(
    ...     chunk_size=3,
    ...     starts=[0, 1, 2],
    ...     stops=[5, 4, 3],
    ...     confidence_matrix=conf_mat)
    >>> len(supp_mat) == len(conf_mat)
    True
    >>> round(supp_mat[1, 1], 4)
    0.7
    >>> round(supp_mat[1, 2], 4)
    0.6333
    >>> round(supp_mat[1, 3], 4)
    0.5
    >>> round(supp_mat[1, 4], 4)
    0.5
    >>> (1, 0) in supp_mat
    False
    >>> round(supp_mat[2, 2], 4)
    0.2
    >>> round(supp_mat[2, 3], 4)
    0.2
    >>> (2, 4) in supp_mat
    False
    >>>
    >>> conf_mat = {
    ...     (0, 0): 1.0, (0, 1): 0.1, (0, 2): 0.05, (0, 3): 0.01, (0, 4): 0.01, (0, 5): 0.01, (0, 6): 0.9, (0, 7): 1.0,
    ...     (1, 1): 0.9, (1, 2): 0.75, (1, 3): 0.5, (1, 4): 0.39, (1, 5): 0.1, (1, 6): 0.1,
    ...     (2, 2): 0.2, (2, 3): 0.49, (2, 4): 0.6, (2, 5): 0.89,
    ... }
    >>> supp_mat = fill_support_matrix(
    ...     chunk_size=3,
    ...     starts=[9, 10, 11],
    ...     stops=[16, 15, 14],
    ...     confidence_matrix=conf_mat)
    >>> len(supp_mat) == len(conf_mat)
    True
    >>> supp_mat[1, 1]
    0.825
    >>> round(supp_mat[1, 2], 4)
    0.7167
    >>> round(supp_mat[1, 3], 4)
    0.5467
    >>> round(supp_mat[1, 4], 4)
    0.33
    >>> round(supp_mat[1, 5], 4)
    0.1967
    >>> round(supp_mat[1, 6], 4)
    0.1
    >>> (2, 1) in supp_mat
    False
    >>> round(supp_mat[2, 2], 4)
    0.345
    >>> round(supp_mat[2, 3], 4)
    0.43
    >>> round(supp_mat[2, 4], 4)
    0.66
    >>> round(supp_mat[2, 5], 4)
    0.745
    >>> (2, 6) in supp_mat
    False
    """

    support_matrix: Dict[Tuple[int, int], float] = {}

    half_chunk: int = int((chunk_size - 1) / 2)
    row_count = len(starts)

    # This is the amount by which we need to slide all the start and stop values
    # over to the left to match the column values in the confidence matrix for each
    # subfamily alignment.
    # See the documentation for `start` on the `Alignment` class.
    # FIXME: In some cases this may not be set correctly - how so??
    # This is a problem, the tests don't cover this and they only pass because
    # the value is zero in both cases.
    position_offset = starts[0]  # start of skip state

    # Note: start and stop values are closed ranges, so when a stop value is
    # used in a range() call, the stop value passed to range needs to be
    # incremented by 1. For example, range(start, stop + 1). This convention is
    # consistent through the algorithm below.

    for row_index in range(0, row_count):
        start: int = starts[row_index] - position_offset
        stop: int = stops[row_index] - position_offset

        prev_chunk_start = start
        prev_chunk_stop = stop

        sum_of_scores = 0.0

        for col_index in range(start, stop + 1):
            chunk_start = max(col_index - half_chunk, start)
            chunk_stop = min(col_index + half_chunk, stop)

            column_count = chunk_stop - chunk_start + 1

            if col_index == start:
                # Compute the entire chunk if we're working on the first column
                # in the row. This bootstraps our more efficient summation
                # later.
                sum_of_scores = 0.0
                for sum_index in range(chunk_start, chunk_stop + 1):
                    sum_of_scores += confidence_matrix[row_index, sum_index]
            else:
                # Subtract the value from the left-most column, but only if it
                # has slipped out of our window. Add the value from the
                # right-most column, but only if it is newly part of our window.
                if chunk_start > prev_chunk_start:
                    sum_of_scores -= confidence_matrix[
                        row_index, prev_chunk_start
                    ]
                if chunk_stop > prev_chunk_stop:
                    sum_of_scores += confidence_matrix[row_index, chunk_stop]

            support_matrix[row_index, col_index] = sum_of_scores / column_count

            prev_chunk_start = chunk_start
            prev_chunk_stop = chunk_stop

    return support_matrix
