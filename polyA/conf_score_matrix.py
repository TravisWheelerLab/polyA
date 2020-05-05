from typing import Dict, List, Tuple

from align_score_matrix import AlignScoreMatrix
from confidence_cm import confidence_cm

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
