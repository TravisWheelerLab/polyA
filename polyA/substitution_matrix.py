import re
from typing import Dict, TextIO

SubMatrix = Dict[str, int]


def load_substitution_matrix(file: TextIO) -> SubMatrix:
    """
    Reads the score matrix from a file and returns a dictionary that
    maps '<char1><char2>' to the score from the input substitution
    matrix.

    For example: 'AA' => 8

    TODO: Write a doctest for this
    """
    sub_matrix: SubMatrix = {}

    # TODO: Add all ambiguity codes
    nucleotide_codes = "AGCTYRWSKMDVHBXN."

    for code1 in nucleotide_codes:
        for code2 in nucleotide_codes:
            sub_matrix[code1 + code2] = 0

    line = next(file)
    line = re.sub(r"^\s+", "", line)
    line = re.sub(r"\s+$", "", line)
    chars = re.split(r"\s+", line)

    count: int = 0
    for line in file:
        line = re.sub(r"^\s+", "", line)
        line = re.sub(r"\s+$", "", line)
        sub_scores = re.split(r"\s+", line)
        for i in range(len(sub_scores)):
            sub_matrix[chars[count] + chars[i]] = int(sub_scores[i])
        count += 1

    return sub_matrix
