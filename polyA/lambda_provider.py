import os
import re
from tempfile import NamedTemporaryFile
from typing import Callable, Dict

LambdaProvider = Callable[[str], float]
"""
A callback that accepts a substitution matrix name
and returns the corresponding lambda value for that
matrix.
"""


class ConstantLambdaProvider:
    """
    A test double that returns a fixed value regardless of the
    value passed when called.

    >>> p: LambdaProvider = ConstantLambdaProvider(1.0)
    >>> p("foo")
    1.0
    """

    _lambda: float

    def __init__(self, value: float):
        self._lambda = value

    def __call__(self, name: str) -> float:
        return self._lambda


class EaselLambdaProvider:
    """
    A provider that executes Easel behind the scenes to recover a
    lambda value for the given matrix.
    """

    _path: str
    _matrices: Dict[str, Dict[str, int]]

    def __init__(self, easel_path: str, matrices: Dict[str, Dict[str, int]]):
        self._path = easel_path
        self._matrices = matrices

    def __call__(self, name: str) -> float:
        if name not in self._matrices:
            raise RuntimeError(f"matrix named '{name}' not found")

        # Create the temporary matrix file
        with NamedTemporaryFile("w", delete=False) as matrix_file:
            temp_matrix_path = matrix_file.name
            matrix = self._matrices[name]
            chars = [k[0] for k in matrix]
            matrix_file.write(" ".join(chars) + "\n")
            for char1 in chars:
                line = " ".join([str(matrix[char1 + char2]) for char2 in chars])
                matrix_file.write(line + "\n")

        # Run Easel
        esl_stream = os.popen(
            self._path + "esl_scorematrix --dna " + temp_matrix_path
        )
        esl_output = esl_stream.read()
        esl_output_list = re.split(r"\n+", esl_output)
        lambda_list = re.split(r"\s+", esl_output_list[1])

        # Clean up by removing the temporary matrix file
        os.remove(temp_matrix_path)

        return float(lambda_list[2])
