import math
from logging import getLogger
from timeit import default_timer as timer
from typing import List, Optional, Tuple

from .constants import DEFAULT_CHANGE_PROB, DEFAULT_CHUNK_SIZE
from .origin_matrix import OriginMatrix
from .prob_matrix import ProbMatrix
from .support_matrix import SupportMatrix


def fill_prob_matrix(
    support_matrix: SupportMatrix,
    columns: Optional[List[int]] = None,
    row_count: Optional[int] = None,
    benchmark: bool = False,
    change_prob: float = DEFAULT_CHANGE_PROB,
    chunk_size: int = DEFAULT_CHUNK_SIZE,
) -> Tuple[ProbMatrix, OriginMatrix]:
    """
    Fills in the probability score matrix from the support matrix. Also fills
    the origin matrix for convenience. Returns a tuple that contains the
    probability matrix and the origin matrix, in that order..

    support_matrix -- a column-major nested list (matrix) of support scores
    columns        -- the list of sequence column indices, treated as
                      left-most edges of the chunk window, that should be
                      included
    row_count      -- the total number of rows that will be processed, will
                      be computed from the other arguments if omitted
    benchmark      -- whether or not to collect and output benchmarking data
                      to the logger
    change_prob    -- probability of changing sequence
    chunk_size     -- the width of the "window" analyzed at one time

    We omit filling the first column because we want its probabilities to be
    zero anyway.

    The basic algorithm is described below. All calculations happen in log
    space.

    look at all i's in j-1
        mult by confidence in current cell
        if comes from same i, mult by higher prob
        else - mult by lower prob /(numseqs-1) -> so sum of all probs == 1
     return max
     
    >>> supp_mat = {(0, 0): 0.7, (0, 1): 0.1, (0, 2): 0.1, (1, 0): 0.1, (1, 1): 0.9, (1, 2): 0.9}
    >>> non_cols = [0, 1, 2]
    >>> (prob_mat, orig_mat) = fill_prob_matrix(supp_mat, non_cols, 2)
    >>> orig_mat
    {(0, 0): -1, (1, 0): -1, (0, 1): 0, (1, 1): 1, (0, 2): 0, (1, 2): 1}
    >>> f"{prob_mat[0,0]:.1f}"
    '0.0'
    >>> f"{prob_mat[1,0]:.1f}"
    '0.0'
    >>> f"{prob_mat[0,1]:.4f}"
    '-2.3026'
    >>> f"{prob_mat[1,1]:.4f}"
    '-0.1054'
    >>> f"{prob_mat[0,2]:.4f}"
    '-4.6052'
    >>> f"{prob_mat[1,2]:.4f}"
    '-0.2107'
    
    unrealistic numbers, but want to test when it switches subfams
    
    >>> supp_mat2 = {(0, 0): 10.0, (0, 1): (10 ** -50), (0, 2): (10 ** -50), (1, 0): (10 ** -50), (1, 1): 10.0, (1, 2): 10.0}
    >>> non_cols = [0, 1, 2]
    >>> (prob_mat2, orig_mat2) = fill_prob_matrix(supp_mat2, non_cols, 2)
    >>> orig_mat2
    {(0, 0): -1, (1, 0): -1, (0, 1): 0, (1, 1): 1, (0, 2): 1, (1, 2): 1}
    >>> f"{prob_mat2[0,0]:.1f}"
    '0.0'
    >>> f"{prob_mat2[1,0]:.1f}"
    '0.0'
    >>> f"{prob_mat2[0,1]:.4f}"
    '-115.1293'
    >>> f"{prob_mat2[1,1]:.4f}"
    '2.3026'
    >>> f"{prob_mat2[0,2]:.4f}"
    '-216.4430'
    >>> f"{prob_mat2[1,2]:.4f}"
    '4.6052'
	
    """
    same_prob = 1.0 - change_prob

    logged_change_prob = math.log(change_prob)
    logged_same_prob = math.log(same_prob)

    if benchmark:
        start = timer()

    origins: OriginMatrix = {}
    probabilities: ProbMatrix = {}

    if columns is None or row_count is None:
        maxCol = 0
        maxRow = 0
        for (row, col) in support_matrix:
            if col > maxCol:
                maxCol = col
            if row > maxRow:
                maxRow = row

    if columns is None:
        columns = list(range(maxCol + 1))

    if row_count is None:
        row_count = maxRow + 1

    for rowIndex in range(row_count):
        origins[(rowIndex, 0)] = -1
        probabilities[(rowIndex, 0)] = 0.0

    # colIndex is j in the Perl version
    colIndex = 1

    # Skip the first column because we already filled it above.
    # We need to run off the end of the last column by `chunk_size` since we
    # slide the window along one column at a time. So we augment the selected
    # columns with the necessary additional columns.
    # colIndexIndex is col in the Perl version
    for colIndexIndex in range(1, len(columns) + chunk_size - 1):
        inProvidedColumns = colIndexIndex < len(columns)

        # Matches behavior on lines 819-823
        if inProvidedColumns:
            colIndex = columns[colIndexIndex]
        else:
            colIndex += 1

        for rowIndex in range(row_count):
            maxScore = -math.inf
            maxScoreIndex = 0

            # Mapped from line 831
            support = support_matrix.get((rowIndex, colIndex), None)
            if support is not None and support != 0:
                supportLog = math.log(support)

                for innerRowIndex in range(row_count):
                    # Omitting check on line 841 because we expanded columns
                    if inProvidedColumns:
                        probValue = probabilities.get(
                            (innerRowIndex, columns[colIndexIndex - 1]), None
                        )
                    else:
                        probValue = probabilities.get(
                            (innerRowIndex, colIndex - 1), None
                        )

                    if probValue is not None:
                        score = supportLog + probValue

                        if rowIndex == innerRowIndex:
                            # Same sequence (higher probability)
                            score += logged_same_prob
                        else:
                            # Different sequence (lower probability)
                            score += logged_change_prob

                        # Update our current max score
                        if score > maxScore:
                            maxScore = score
                            maxScoreIndex = innerRowIndex

                probabilities[(rowIndex, colIndex)] = maxScore
                origins[(rowIndex, colIndex)] = maxScoreIndex

    if benchmark:
        end = timer()
        logger = getLogger(__name__)
        logger.info(f"benchmark: {end - start}s")

    return probabilities, origins
