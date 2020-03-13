from typing import Iterable, List
from .align_score_matrix import AlignScoreMatrix


# TODO: If we improved the AlignScoreMatrix type we could get column and row counts from that
def non_empty_columns(
    align_matrix: AlignScoreMatrix, seq_length: int, seq_count: int
) -> Iterable[int]:
    """
    TODO: Better docstring
    # puts all columns that are not empty into @Columns, so when
    # I loop through hash I can use the
    # vals in @Columns - this will skip over empty columns

    >>> a = {(0, 0): 1.0, (0, 1): 1.0}
    >>> non_empty_columns(a, 3, 1)
    [0, 1]
    """
    columns: List[int] = []
    for j in range(seq_length):
        is_empty = True
        for i in range(seq_count):
            if (i, j) in align_matrix:
                is_empty = False
                break

        if not is_empty:
            columns.append(j)

    return columns
