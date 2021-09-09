from typing import Dict, List
import re


def confidence_cm(
    region: List[float],
    subfam_counts: Dict[str, float],
    subfams: List[str],
    subfam_rows: List[int],
    repeats: int,
    node_confidence_bool: bool,
) -> List[float]:
    """
    Computes confidence values for competing annotations using alignment and tandem
    repeat scores. Loops through the array once to find sum of 2^every_hit_score in
    region, then loops back through to calculate confidence. Converts the alignment
    score to account for lambda before summing.

    If command line option for subfam_counts, this is included in confidence math.

    input:
    region: list of scores for competing annotations
    subfam_counts: dict that maps subfam name to it's count info
    subfams: list of subfam names
    subfam_rows: list of subfamily rows that correspond to the the subfams of the region scores
    repeats: number of tandem repeat scores found at the end of the region
    node_confidence_bool: if False (0) - confidence for filling DP matrices, if True - confidence for nodes

    output:
    confidence_list: list of confidence values for competing annotations, each input alignment
    and tandem repeat score will have one output confidence score

    >>> counts = {"s1": .33, "s2": .33, "s3": .33}
    >>> subs = ["s1", "s2", "s3"]
    >>> conf = confidence_cm([2, 1, 1], counts, subs, [0, 1, 2], 0, False)
    >>> f"{conf[0]:.2f}"
    '0.50'
    >>> f"{conf[1]:.2f}"
    '0.25'
    >>> f"{conf[2]:.2f}"
    '0.25'

    >>> conf = confidence_cm([0, 100, 100], 0, subs, [0, 1, 2], 0, False)
    >>> f"{conf[0]:.2f}"
    '0.01'

    >>> conf = confidence_cm([0, 100, 100], 0, subs, [0, 1, 2], 0, True)
    >>> f"{conf[0]:.2f}"
    '0.00'

    >>> counts = {"s1": .31, "s2": .31, "s3": .31, "Tandem Repeat": .06}
    >>> subs = ["s1", "s2", "s3", "Tandem Repeat"]
    >>> conf = confidence_cm([2, 1, 0.7], counts, subs, [0, 1, 3], 1, False)
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
    if subfam_counts:
        # alignment scores
        for index in range(len(region) - repeats):
            subfam: str = subfams[subfam_rows[index]]
            m = re.search(r"(.+?)#.+", subfams[subfam_rows[index]])
            if m:
                subfam = m.group(1)
            converted_score = (2 ** int(region[index])) * subfam_counts[subfam]
            confidence_list.append(converted_score)
            score_total += converted_score
        # TR scores
        for index in range(len(region) - repeats, len(region)):
            subfam = subfams[subfam_rows[index]]
            m = re.search(r"(.+?)#.+", subfams[subfam_rows[index]])
            if m:
                subfam = m.group(1)
            tr_score = (2 ** int(region[index])) * subfam_counts[subfam]
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

    # if skip state confidence is < 1 %, increase it to 1 % and normalize all
    # others do not do this when computing node confidence (skip state is not
    # used)
    if confidence_list[0] < 0.01 and not node_confidence_bool:
        summ = 0.0
        for i in range(1, len(confidence_list)):
            summ += confidence_list[i]

        confidence_list[0] = 0.01
        for i in range(1, len(confidence_list)):
            confidence_list[i] = confidence_list[i] * 0.99 / summ

    return confidence_list


def confidence_only(
    region: List[int],
    lambs: List[float],
) -> List[float]:
    """
    Computes confidence values for competing annotations using alignment scores. Loops
    through the array once to find sum of 2^every_hit_score in
    region, then loops back through to calculate confidence. Converts the alignment
    score to account for lambda before summing.

    input:
    lambdaa: lambda for score matrix used
    region: list of scores for competing annotations
    subfams: list of subfam names

    output:
    confidence_list: list of confidence values for competing annotations

    >>> reg = [100, 55, 1]
    >>> lambs = [.1227] * 3
    >>> conf = confidence_only(reg, lambs)
    >>> conf
    [0.9843787551069454, 0.015380918048546022, 0.0002403268445085316]
    """

    confidence_list: List[float] = []
    score_total: int = 0

    # alignment scores
    for index in range(len(region)):
        converted_score = 2 ** int(region[index] * lambs[index])
        confidence_list.append(converted_score)
        score_total += converted_score

    for index in range(len(region)):
        confidence_list[index] = confidence_list[index] / score_total

    return confidence_list
