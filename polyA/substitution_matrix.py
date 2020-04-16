from typing import Dict, List, TextIO, Tuple

SubstitutionMatrix = Dict[Tuple[str, str], float]
"""
We're using a tuple for the keys so that we don't have to
do anything weird to get a value out of the matrix based on
two characters.

TODO: Explain what this is...
"""


def load_substitution_matrix(
    file: TextIO,
) -> SubstitutionMatrix:
    """
    Load the substitution matrix along with the map of character positions.

    TODO: Explain this better, like what it's used for
    TODO: Make the doctest use realistic input values

    >>> from io import StringIO
    >>> file = StringIO("a b\\n1.0 2.0\\n3.0 4.0")
    >>> s = load_substitution_matrix(file)
    >>> s
    {('a', 'a'): 1.0, ('a', 'b'): 2.0, ('b', 'a'): 3.0, ('b', 'b'): 4.0}
    """
    first_line: str = next(file).replace(" ", "").replace("\t", "").strip()
    labels: List[str] = list(first_line)
    substitution_matrix: SubstitutionMatrix = {}
    for row_index, line in enumerate(file):
        row_scores = line.split(" ")
        for col_index, score in enumerate(row_scores):
            row_label = labels[row_index]
            col_label = labels[col_index]
            substitution_matrix[row_label, col_label] = float(score)

    return substitution_matrix


if __name__ == "__main__":
    import doctest
    doctest.testmod()
