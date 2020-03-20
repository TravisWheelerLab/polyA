from polyA import load_substitution_matrix
from pytest import approx, mark
from typing import List, Tuple


@mark.sub_matrix
def test_load_substitution_matrix():
    file = open("fixtures/sub-matrix.in", "r")
    charPositions, subMatrix = load_substitution_matrix(file)

    assert charPositions["a"] == 0
    assert charPositions["b"] == 1
    assert charPositions["c"] == 2

    assert subMatrix[0] == 1.0
    assert subMatrix[1] == 2.0
    assert subMatrix[2] == 3.0
    assert subMatrix[3] == 4.0
    assert subMatrix[4] == 5.0
    assert subMatrix[5] == 6.0
