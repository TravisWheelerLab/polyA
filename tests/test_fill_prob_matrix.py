from functools import lru_cache
from os import listdir
from os.path import isfile
from polyA import (
    OriginMatrix,
    ProbMatrix,
    deserialize_prob_matrix,
    deserialize_support_matrix,
    fill_prob_matrix,
)
from pytest import approx, mark
from typing import List, Tuple


# +------------------+
# | Helper Functions |
# +------------------+


class Example:
    actual_prob: ProbMatrix
    expected_prob: ProbMatrix
    actual_origin: OriginMatrix
    expected_origin: OriginMatrix


def compare_origin_matrices(example: Example) -> None:
    # for position in actual:
    #     assert position in expected
    #     assert expected[position] == actual[position]
    for position in example.expected_origin:
        assert position in example.actual_origin
        assert example.expected_origin[position] == example.actual_origin[position]


def compare_prob_matrices(example: Example) -> None:
    # for position in actual:
    #     assert position in expected
    #     assert expected[position] == approx(actual[position])
    for position in example.expected_prob:
        assert position in example.actual_prob
        assert example.expected_prob[position] == approx(
            example.actual_prob[position], rel=0.1
        )


@lru_cache(100)
def load_example(id: int) -> Example:
    example = Example()

    with open(f"fixtures/ex{id}.in", "r") as inFile:
        supportLines = inFile.readlines()
        support = deserialize_support_matrix(supportLines)
    example.actual_prob, example.actual_origin = fill_prob_matrix(support)

    with open(f"fixtures/ex{id}.out", "r") as outFile:
        expectedLines = outFile.readlines()
        example.expected_prob, example.expected_origin = deserialize_prob_matrix(
            expectedLines
        )

    return example


# +------------+
# | Test Cases |
# +------------+


# This is the number of "exN.in" and "exN.out" file pairs we have in the
# fixtures/ directory. This is kind of a stupid thing to do, but it is better
# than having to keep track manually.
EXAMPLE_COUNT = 0
for name in listdir("fixtures/"):
    if name.startswith("ex") and name.endswith(".in"):
        EXAMPLE_COUNT += 1


@mark.slow
@mark.prob_matrix
@mark.parametrize("index", range(EXAMPLE_COUNT))
def test_origin_example(index: int):
    example = load_example(index)
    compare_origin_matrices(example)


@mark.slow
@mark.origin_matrix
@mark.parametrize("index", range(EXAMPLE_COUNT))
def test_prob_example(index: int):
    example = load_example(index)
    compare_prob_matrices(example)
