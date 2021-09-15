from typing import Dict


def sub_repeat_scores(
    start: int, end: int, repeat_scores: Dict[int, float]
) -> float:
    """
    Sums tandem repeat scores from columns start to end.

    input:
    start: column index to start sum
    end: column index to end sum
    repeat_scores: hash implementation of 1d array that maps column index
    to tandem repeat score

    output:
    summ: sum of the tandem repeat scores
    """

    summ: float = 0.0
    curr_index = start
    while curr_index < end:
        if curr_index in repeat_scores:
            summ += repeat_scores[curr_index]
        curr_index += 1
    return summ
