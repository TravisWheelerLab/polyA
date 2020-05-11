from typing import Dict, List, TextIO, Tuple

SubstitutionMatrix = Dict[Tuple[str, str], float]
"""
We're using a tuple for the keys so that we don't have to
do anything weird to get a value out of the matrix based on
two characters.

Score matrix used when calculating alignment scores. Holds values to 
give when aligning 2 letters (nucleotides). 

Ex: VERY basic substitution matrix will give a score of 1 for matchinn nucleotides and 
0 for mismatching nucleotides
	A	C	T	G
A	1	0	0	0
C	0	1	0	0
T	0 	0	1	0
G	0	0	0	1

"""


def load_substitution_matrix(file: TextIO,) -> SubstitutionMatrix:
    """
    Load the substitution matrix along with the map of character positions.

    User will input substitution matrix that is used for scoring alignments.
    If lambda is not given, it will be calculated based on this substitution matrix.
    
    >>> from io import StringIO
    >>> file = StringIO("A C\\n1 0\\n0 1\\n")
    >>> s = load_substitution_matrix(file)
    >>> s
    {('A', 'A'): 1, ('A', 'C'): 0, ('C', 'A'): 0, ('C', 'C'): 1}
    """
    first_line: str = next(file).replace(" ", "").replace("\t", "").strip()
    labels: List[str] = list(first_line)
    substitution_matrix: SubstitutionMatrix = {}
    for row_index, line in enumerate(file):
        row_scores = line.split(" ")
        for col_index, score in enumerate(row_scores):
            row_label = labels[row_index]
            col_label = labels[col_index]
            substitution_matrix[row_label, col_label] = int(score)

    return substitution_matrix


if __name__ == "__main__":
    import doctest

    doctest.testmod()
