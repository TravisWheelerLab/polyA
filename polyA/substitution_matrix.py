import re
import logging
from typing import Dict, List, Optional, TextIO, Tuple

from .lambda_provider import LambdaProvider, ConstantLambdaProvider


_logger = logging.root.getChild(__name__)


class ParseError(RuntimeError):
    pass


ScoresDict = Dict[str, int]


class SubMatrix:
    """
    A single substitution matrix, including its name, lambda value, and assumed
    character background frequencies.
    The lambda value may be `None`, in which case it must be computed
    before the matrix can be used.
    The background frequencies may also be 'None' if complexity adjusted scoring
    will not be used.

    A single matrix maps '<char1><char2>' (a pair of characters from
    the input alphabet) to the score from the input substitution
    matrix.

    For example: 'AA' => 8
    """

    lamb: float
    name: str
    scores: ScoresDict
    background_freqs: Optional[Dict[str, float]]

    def __init__(self, name: str, lamb: float):
        self.lamb = lamb
        self.name = name
        self.scores = {}
        self.background_freqs = None


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


def _parse_background_freqs(line: str) -> Dict[str, float]:
    """
    Attempt to parse a matrix background frequencies from the given line, containing
    each char and it's assumed background frequency.
    """
    matrix_background_freqs: Dict[str, float] = {}
    if line.strip().upper().startswith("BACKGROUND FREQS"):
        clean_line = re.sub(r"^\s+|\s+$", "", line)
        if line.strip().upper().startswith("BACKGROUND FREQS"):
            # convert to dictionary
            string_dict = clean_line[
                clean_line.find("{") : clean_line.find("}") + 1
            ]
            matrix_background_freqs = eval(string_dict)
    return matrix_background_freqs


def _parse_chars(line: str) -> List[str]:
    clean_line = re.sub(r"^\s+|\s+$", "", line)
    return re.split(r"\s+", clean_line)


def load_substitution_matrices(
    file: TextIO,
    lambda_provider: LambdaProvider,
    complexity_adjustment: bool,
    alphabet_chars: str = "AGCTYRWSKMDVHBXN",
    ambiguity_chars: str = "-.",
) -> SubMatrixCollection:
    """
    Reads a set of score matrices from a file and returns a
    `SubMatrixCollection` that contains the matrices and any lambda
    values associated with them.

    Each matrix has a name, which is then specified for each
    alignment that is to be adjudicated.

    >>> _lamb_provider = ConstantLambdaProvider(0.1227)
    >>> matrices_file = "fixtures/ultra_test_files/ex13.fa.cm.matrix"
    >>> with open(matrices_file) as _sub_matrices_file:
    ...     sub_matrices = load_substitution_matrices(_sub_matrices_file, _lamb_provider, False)
    >>> len(sub_matrices)
    1
    >>> "matrix1" in sub_matrices
    True
    >>> matrix = sub_matrices["matrix1"]
    >>> matrix.lamb
    0.1227
    >>> matrix.name
    'matrix1'
    >>> matrix.scores['AA']
    8
    >>> matrix.scores['AT']
    -15
    >>> matrix.scores['TT']
    8
    >>> matrix.scores['GA']
    -2
    >>> matrix.scores['GC']
    -13
    >>> matrix.scores['CC']
    10
    >>> matrix.background_freqs is None
    True
    >>> with open(matrices_file) as _sub_matrices_file:
    ...     sub_matrices = load_substitution_matrices(_sub_matrices_file, _lamb_provider, True)
    >>> matrix = sub_matrices["matrix1"]
    >>> matrix.background_freqs is None
    False
    >>> matrix.background_freqs
    {'A': 0.295, 'G': 0.205, 'C': 0.205, 'T': 0.295}
    """
    _logger.debug(f"load_substitution_matrix({file.name})")

    collection: SubMatrixCollection = {}

    while True:
        next_scores: ScoresDict = {}
        matrix_background_freqs: Dict[str, float] = {}

        for char1 in alphabet_chars + ambiguity_chars:
            for char2 in alphabet_chars + ambiguity_chars:
                next_scores[char1 + char2] = 0

        try:
            header_line = next(file)
            name, lamb = _parse_matrix_header(header_line)
        except StopIteration:
            break

        try:
            unknown_line = next(
                file
            )  # this could be background freqs or matrix chars
            if unknown_line.strip().upper().startswith("BACKGROUND FREQS"):
                if complexity_adjustment:
                    matrix_background_freqs = _parse_background_freqs(
                        unknown_line
                    )
                    if len(matrix_background_freqs) == 0:
                        # no matrix background frequencies were found
                        raise ParseError(
                            f"cannot use complexity adjusted scoring, missing background frequencies: '{file.name}'"
                        )
                try:
                    # skip background freqs line
                    # this should be the substitution matrix chars
                    unknown_line = next(file)
                except StopIteration:
                    raise ParseError("dangling substitution matrix header")
        except StopIteration:
            break

        chars_line = unknown_line
        chars = _parse_chars(chars_line)

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
        if complexity_adjustment:
            sub_matrix.background_freqs = matrix_background_freqs
    return collection
