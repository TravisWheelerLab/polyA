from argparse import ArgumentError
import os
import re
import subprocess
from tempfile import NamedTemporaryFile
from typing import Callable, Dict, Optional

from .performance import timeit

LambdaProvider = Callable[[Dict[str, int]], float]
"""
A callback that accepts a substitution matrix name
and returns the corresponding lambda value for that
matrix.
"""


class LambdaProviderException(Exception):
    return_code: int
    stdout: str
    stderr: str

    def __init__(self, return_code: int, stdout: str, stderr: str):
        self.return_code = return_code
        self.stdout = stdout
        self.stderr = stderr


class ConstantLambdaProvider:
    """
    A test double that returns a fixed value regardless of the
    value passed when called.

    >>> p: LambdaProvider = ConstantLambdaProvider(1.0)
    >>> p({})
    1.0
    """

    _lambda: float

    def __init__(self, value: float):
        self._lambda = value

    def __call__(self, matrix: Dict[str, int]) -> float:
        return self._lambda


class EaselLambdaProvider:
    """
    A provider that executes Easel behind the scenes to recover a
    lambda value for the given matrix.
    """

    _path: str

    def __init__(self, easel_path: str):
        self._path = easel_path

    @timeit()
    def __call__(self, matrix: Dict[str, int]) -> float:
        # Create the temporary matrix file
        with NamedTemporaryFile("w", delete=False) as matrix_file:
            temp_matrix_path = matrix_file.name
            chars = set([k[0] for k in matrix])
            matrix_file.write(" ".join(chars) + "\n")

            for char1 in chars:
                line = " ".join([str(matrix[char1 + char2]) for char2 in chars])
                matrix_file.write(line + "\n")

        # Run Easel
        esl_process = subprocess.run(
            [self._path, "--dna", temp_matrix_path],
            capture_output=True,
            text=True,
        )

        if esl_process.returncode != 0:
            raise LambdaProviderException(
                esl_process.returncode, esl_process.stdout, esl_process.stderr
            )

        esl_output = esl_process.stdout
        esl_output_list = re.split(r"\n+", esl_output)
        lambda_list = re.split(r"\s+", esl_output_list[1])

        # Clean up by removing the temporary matrix file
        os.remove(temp_matrix_path)

        return float(lambda_list[2])
