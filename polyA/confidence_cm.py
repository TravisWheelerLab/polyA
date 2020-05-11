from typing import List


def confidence_cm(lambda_value: float, region: List[float]) -> str:
    """    
    input: list of scores for competing annotations, lambda value for score matrix used
    output: list of confidence values for competing annotations, each input alignment 
    score will have one output confidence score
    
    computes confidence values for competing annoations using alignment scores
    Loops through the array once to find sum, then loops back through to calculate confidence

    TODO: Can this return a list instead of a string?  -- YES (from Kaitlin)

    >>> confidence_cm(0.5, [1.0, 0.0, -1.0])
    '1.0 0.0 0.0'
    """
    converted_scores = [
        2 ** (score / lambda_value) if score > 0 else 0.0 for score in region
    ]
    score_total = sum(converted_scores)
    confidence_values = (str(score / score_total) for score in converted_scores)
    return " ".join(confidence_values)


if __name__ == "__main__":
    import doctest

    doctest.testmod()
