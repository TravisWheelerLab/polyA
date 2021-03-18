import re
import logging
from typing import Dict, List, Optional, TextIO, Tuple

from .lambda_provider import LambdaProvider


_logger = logging.root.getChild(__name__)


class ParseError(RuntimeError):
    pass


ScoresDict = Dict[str, int]


class SubMatrix:
    """
    A single substitution matrix, including its name and lambda value.
    The lambda value may be `None`, in which case it must be computed
    before the matrix can be used.

    A single matrix maps '<char1><char2>' (a pair of characters from
    the input alphabet) to the score from the input substitution
    matrix.

    For example: 'AA' => 8
    """

    lamb: float
    name: str
    scores: ScoresDict

    def __init__(self, name: str, lamb: float):
        self.lamb = lamb
        self.name = name
        self.scores = {}


SubMatrixCollection = Dict[str, SubMatrix]


def _parse_matrix_header(line: str) -> Tuple[str, Optional[float]]:
    """
    Attempt to parse a matrix header from the given line, containing
    a matrix name and, optionally, a lambda value.
    """
    clean_line = re.sub(r"^\s+|\s+$", "", line)
    tokens = re.split(r"\s", clean_line)

    if len(tokens) == 0:
        raise ParseError(f"incomplete substitution matrix header: '{line}'")
    if len(tokens) > 2:
        raise ParseError(f"invalid substitution matrix header: '{line}'")

    if len(tokens) == 1:
        return tokens[0], None
    # len(tokens) == 2
    return tokens[0], float(tokens[1])


def _parse_chars(line: str) -> List[str]:
    clean_line = re.sub(r"^\s+|\s+$", "", line)
    return re.split(r"\s+", clean_line)


def load_substitution_matrices(
    file: TextIO,
    lambda_provider: LambdaProvider,
    alphabet_chars: str = "AGCTYRWSKMDVHBXN",
    ambiguity_chars: str = ".",
) -> SubMatrixCollection:
    """
    Reads a set of score matrices from a file and returns a
    `SubMatrixCollection` that contains the matrices and any lambda
    values associated with them.

    Each matrix has a name, which is then specified for each
    alignment that is to be adjudicated.

    TODO: Write a doctest for this
    TODO: Add additional ambiguity codes to default value
    """
    _logger.debug(f"load_substitution_matrix({file.name})")

    collection: SubMatrixCollection = {}

    while True:
        next_scores: ScoresDict = {}

        for char1 in alphabet_chars + ambiguity_chars:
            for char2 in alphabet_chars + ambiguity_chars:
                next_scores[char1 + char2] = 0

        try:
            header_line = next(file)
            name, lamb = _parse_matrix_header(header_line)
        except StopIteration:
            break

        try:
            chars_line = next(file)
            chars = _parse_chars(chars_line)
        except StopIteration:
            raise ParseError("dangling substitution matrix header")

        count: int = 0
        for line in file:
            if "//" in line:
                break

            clean_line = re.sub(r"^\s+|\s+$", "", line)
            sub_scores = re.split(r"\s+", clean_line)
            for i in range(len(sub_scores)):
                next_scores[chars[count] + chars[i]] = int(sub_scores[i])
            count += 1

        if len(next_scores) == 0:
            raise ParseError(
                f"missing substitution matrix values: '{file.name}'"
            )

        if lamb is None:
            lamb = lambda_provider(next_scores)

        sub_matrix = SubMatrix(name, lamb)
        sub_matrix.scores = next_scores
        collection[name] = sub_matrix

    return collection
