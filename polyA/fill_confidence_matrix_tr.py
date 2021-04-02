from typing import Dict, List, Tuple

from .fill_confidence_matrix import fill_confidence_matrix
from .performance import timeit


@timeit
def fill_confidence_matrix_tr(
    columns: List[int],
    subfam_counts: Dict[str, float],
    subfams: List[str],
    active_cells: Dict[int, List[int]],
    align_matrix: Dict[Tuple[int, int], float],
    repeat_scores: Dict[int, float],
) -> Dict[Tuple[int, int], float]:
    """
    Mostly identical to fill_confidence_matrix, but takes into account that some
    columns in the alignment matrix may hold tandem repeats.

    Inputs:

    Superset of the inputs to fill_confidence_matrix. Additions are listed
    below.

    repeat_scores - maps column indices to tandem repeat scores.

    Outputs:

    Same as `fill_confidence_matrix`.
    """
    confidence_matrix = fill_confidence_matrix(
        columns,
        subfam_counts,
        subfams,
        active_cells,
        align_matrix,
    )
    tr_confidence_matrix = fill_confidence_matrix(
        repeat_scores.keys(),
        subfam_counts,
        subfams,
        active_cells,
        align_matrix,
    )

    confidence_matrix.update(tr_confidence_matrix)

    return confidence_matrix
