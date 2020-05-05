from typing import List


def confidence_cm(lambda_value: float, region: List[float]) -> str:
    """
    Loops through the array once to convert scores then produces
    a string of scores.

    TODO (Kaitlin): Explain this a little better

    >>> confidence_cm(0.5, [1.0, 0.0, -1.0])
    '1.0 0.0 0.0'
    """
    converted_scores = [2 ** (score / lambda_value) if score > 0 else 0.0 for score in region]
    score_total = sum(converted_scores)
    confidence_values = (str(score / score_total) for score in converted_scores)
    return " ".join(confidence_values)


if __name__ == "__main__":
    import doctest

    doctest.testmod()
