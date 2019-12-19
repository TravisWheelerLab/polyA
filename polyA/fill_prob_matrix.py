import math
from typing import Any
from typing import Dict, List, Tuple

from .constants import CHUNK_SIZE, LOGGED_CHANGE_PROB, LOGGED_SAME_PROB


def fill_prob_matrix(
    supports: Dict[Tuple[int, int], float], columns: List[int], rowCount: int
) -> Tuple[Any, Any]:
    """
    Fills in the probability score matrix from the support matrix. Also fills
    the origin matrix for convenience.

    supports -- a column-major nested list (matrix) of support scores
    columns  -- the list of sequence column indices, treated as left-most
                edges of the chunk window, that should be included

    We skip the first column because we want its probabilities to be zero
    anyway.

    The basic algorithm is described below. All calculations happen in log
    space.

    look at all i's in j-1
        mult by confidence in current cell
        if comes from same i, mult by higher prob
        else - mult by lower prob /(numseqs-1) -> so sum of all probs == 1
     return max
    """
    probabilities: Dict[Tuple[int, int], float] = {}
    origins: Dict[Tuple[int, int], int] = {}

    colCount = len(columns) + CHUNK_SIZE - 1
    for colIndex in range(1, colCount):
        # Omitting lines 819-823 here so we can verify their necessity

        for rowIndex in range(rowCount):
            maxScore = -math.inf
            maxScoreIndex = 0

            supportLog = -math.inf  # Not sure what this is for (line 832)

            if True:  # Figure out how we will map these (line 831)
                for innerRowIndex in range(rowCount):
                    score = 0.0  # Initial value doesn't matter
                    exists = False

                    # Omitting check on line 841 to verify its necessity
                    cellValue = probabilities.get(
                        (innerRowIndex, columns[colIndex - 1]), default=None
                    )
                    if cellValue is not None:
                        score = supportLog + cellValue
                        exists = True

                    if exists:
                        if rowIndex == innerRowIndex:
                            # Same sequence (higher probability)
                            score = score + LOGGED_SAME_PROB
                        else:
                            # Different sequence (lower probability)
                            score = score + LOGGED_CHANGE_PROB

                        # Update our current max score
                        if score > maxScore:
                            maxScore = score
                            maxScoreIndex = innerRowIndex

                probabilities[(rowIndex, colIndex)] = maxScore
                origins[(rowIndex, colIndex)] = maxScoreIndex

    return (probabilities, origins)
