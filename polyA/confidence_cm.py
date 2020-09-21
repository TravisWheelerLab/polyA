from typing import Dict, List


def confidence_cm(
    infile: str,
    region: List[float],
    subfam_counts: Dict[str, float],
    subfams: List[str],
    subfam_rows: List[int],
    repeats: int,
) -> List[float]:
    """
    Computes confidence values for competing annotations using alignment and tandem
    repeat scores. Loops through the array once to find sum of 2^every_hit_score in
    region, then loops back through to calculate confidence. Converts the alignment
    score to account for lambda before summing.

    If command line option for subfam_counts, this is included in confidence math.

    input:
    lambdaa: lambda value for input sub_matrix (scaling factor)
    infile: test if subfam_counts infile included at command line
    region: list of scores for competing annotations
    subfam_counts: dict that maps subfam name to it's count info
    subfams: list of subfam names
    subfam_rows: list of subfamily rows that correspond to the the subfams of the region scores
    repeats: number of tandem repeat scores found at the end of the region

    output:
    confidence_list: list of confidence values for competing annotations, each input alignment
    and tandem repeat score will have one output confidence score

    >>> counts = {"s1": .33, "s2": .33, "s3": .33}
    >>> subs = ["s1", "s2", "s3"]
    >>> conf = confidence_cm("infile", [2, 1, 1], counts, subs, [0, 1, 2], 0)
    >>> f"{conf[0]:.2f}"
    '0.50'
    >>> f"{conf[1]:.2f}"
    '0.25'
    >>> f"{conf[2]:.2f}"
    '0.25'

    >>> counts = {"s1": .31, "s2": .31, "s3": .31, "Tandem Repeat": .06}
    >>> subs = ["s1", "s2", "s3", "Tandem Repeat"]
    >>> conf = confidence_cm("infile", [2, 1, 0.7], counts, subs, [0, 1, 3], 1)
    >>> f"{conf[0]:.2f}"
    '0.65'
    >>> f"{conf[1]:.2f}"
    '0.32'
    >>> f"{conf[2]:.2f}"
    '0.03'
    """

    confidence_list: List[float] = []
    score_total: int = 0

    # if command line option to include subfam_counts
    if infile:
        # alignment scores
        for index in range(len(region) - repeats):
            converted_score = (
                2 ** int(region[index])
            ) * subfam_counts[subfams[subfam_rows[index]]]
            confidence_list.append(converted_score)
            score_total += converted_score
        # TR scores
        for index in range(len(region) - repeats, len(region)):
            tr_score = (2 ** int(region[index])) * subfam_counts[
                subfams[subfam_rows[index]]
            ]
            confidence_list.append(tr_score)
            score_total += tr_score

    # don't include subfam counts (default)
    else:
        # alignment scores
        for index in range(len(region) - repeats):
            converted_score = 2 ** int(region[index])
            confidence_list.append(converted_score)
            score_total += converted_score
        # TR scores
        for index in range(len(region) - repeats, len(region)):
            tr_score = 2 ** int(region[index])
            confidence_list.append(tr_score)
            score_total += tr_score

    for index in range(len(region)):
        confidence_list[index] = confidence_list[index] / score_total

    return confidence_list
