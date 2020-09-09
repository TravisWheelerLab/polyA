from typing import Dict, List


# TODO: Rename to be more descriptive
# TODO: Change the type of `infile` to make more sense
def confidence_cm(
    lambdaa: float,
    infile: str,
    region: List[float],
    subfam_counts: Dict[str, float],
    subfams: List[str],
) -> List[float]:
    """
    computes confidence values for competing annotations using alignment scores.
    Loops through the array once to find sum of 2^every_hit_score in region, then
    loops back through to calculate confidence. Converts the score to account for
    lambda before summing.

    If command line option for subfam_counts, this is included in confidence math.

    input:
    lambdaa: lambda value for input sub_matrix (scaling factor)
    infile: test if subfam_counts infile included at command line
    region: list of scores for competing annotations
    subfam_counts: dict that maps subfam name to it's count info
    subfams: list of subfam names

    output:
    confidence_list: list of confidence values for competing annotations, each input alignment
    score will have one output confidence score

    >>> counts = {"s1": .33, "s2": .33, "s3": .33}
    >>> subs = ["s1", "s2", "s3"]
    >>> conf = confidence_cm(0.5, "infile", [2, 1, 1], counts, subs)
    >>> f"{conf[0]:.2f}"
    '0.41'
    >>> f"{conf[1]:.2f}"
    '0.29'
    >>> f"{conf[2]:.2f}"
    '0.29'
    """
    confidence_list: List[float] = []

    score_total: int = 0

    # if command line option to include subfam_counts
    if infile: #FIXME - doesn't always get the correct subfam_counts[] entry because this uses index of the smooshed array instead of the regular one
        for index in range(len(region)):
            converted_score = (2 ** int(region[index] * lambdaa)) * subfam_counts[
                subfams[index]
            ]
            confidence_list.append(converted_score)
            score_total += converted_score

        for index in range(len(region)):
            confidence_list[index] = confidence_list[index] / score_total

    # don't include subfam counts (default)
    else:
        for index in range(len(region)):
            converted_score = 2 ** int(region[index] * lambdaa)
            confidence_list.append(converted_score)
            score_total += converted_score

        for index in range(len(region)):
            confidence_list[index] = confidence_list[index] / score_total

    return confidence_list
