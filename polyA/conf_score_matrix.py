from typing import Dict, List, Tuple

from .align_score_matrix import AlignScoreMatrix
from .confidence_cm import confidence_cm

ConfScoreMatrix = Dict[Tuple[int, int], float]
"""
Typedef to represent a sparse confidence matrix implemented as a
dictionary (for now).
"""


def fill_conf_score_matrix(
    align_score_matrix: AlignScoreMatrix,
    non_empty_columns: List[int],
    n_rows: int,
    lambda_value: float,
) -> ConfScoreMatrix:

    """
    fills confidence matrix from alignment matrix.
    each column in the alignment matrix is a group of competing annotations that are
    input into confidence_cm, the output confidence values are used to populate the confidence 
    matrix    
    
    >>> align_mat = {(0, 0): 0, (0, 1): 100, (0, 2): 99, (1, 0): 100, (1, 1): 100, (1, 2): 100}
    >>> non_cols = [0, 1, 2]
    >>> conf_mat = {}
    >>> conf_mat = fill_conf_score_matrix(align_mat, non_cols, 2, 0.1227)
    >>> print(f"{conf_mat[0,0]:.4f}")
    0.0000
    >>> print(f"{conf_mat[1,0]:.4f}")
    1.0000
    >>> print(f"{conf_mat[0,1]:.4f}")
    0.5000
    >>> print(f"{conf_mat[1,1]:.4f}")
    0.5000
    >>> print(f"{conf_mat[0,2]:.4f}")
    0.0035
    >>> print(f"{conf_mat[1,2]:.4f}")
    0.9965
    
    """
    
    conf_score_matrix: ConfScoreMatrix = {}

    for col in non_empty_columns:
        scores: List[float] = []

        for row in range(n_rows):
            if (row, col) in align_score_matrix:
                scores.append(align_score_matrix[row, col])
            else:
                scores.append(0.0)

        # TODO: Can we just return a list from confidence_cm?
        confidence_temp = confidence_cm(lambda_value, scores).split(" ")

        for row in range(n_rows):
            if confidence_temp[row] != "0":
                conf_score_matrix[row, col] = float(confidence_temp[row])

    return conf_score_matrix
