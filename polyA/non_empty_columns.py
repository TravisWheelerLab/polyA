from typing import Iterable, List
from .align_score_matrix import AlignScoreMatrix


# TODO: If we improved the AlignScoreMatrix type we could get column and row counts from that
def non_empty_columns(
    align_score_matrix: AlignScoreMatrix, seq_length: int, seq_count: int
) -> Iterable[int]:
    """
    Returns an array that holds all non empty columns in AlignScoreMatrix.

    Input alignment file may have empty nucleotide positions (columns) for all alignments.
    The empty columns are not included in non_empty_columns. Instead of looping through
    all columns and jumping over empty ones, we loop through non_empty_columns.
    
    Also used when removing nodes during stitching process. For any node that is
    removed, all corresponding columns are removed from non_empty_columns. When looping
    over non_empty_columns these nodes are ignored when stitching surrounding nodes.

    >>> a = {(0, 0): 1, (0, 1): 1}
    >>> non_empty_columns(a, 3, 1)
    [0, 1]
    """
    columns: List[int] = []
    for j in range(seq_length):
        is_empty = True
        for i in range(seq_count):
            if (i, j) in align_score_matrix:
                is_empty = False
                break

        if not is_empty:
            columns.append(j)

    return columns
