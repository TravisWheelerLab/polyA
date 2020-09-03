import os
import re
from typing import Callable

LambdaProvider = Callable[[], float]


class ConstantLambdaProvider:
    _lambda: float

    def __init__(self, value: float):
        """
        Creates a fake lambda provider that always returns the
        given value.
        """
        self._lambda = value

    def __call__(self) -> float:
        return self._lambda


class EaselLambdaProvider:
    _path: str

    def __init__(self, path: str):
        """
        Creates a new provider instance that runs Easel, which must
        be available at the given path, to produce a value for lambda.
        """
        self._path = path

    def __call__(self) -> float:
        # TODO: How do we know which matrix file to use?
        # Should the matrix file be passed to the call or
        # to the object constructor?

        esl_stream = os.popen(
            self._path + "esl_scorematrix --dna 25p41g_edited.matrix"
        )
        esl_output = esl_stream.read()
        esl_output_list = re.split(r"\n+", esl_output)
        lambda_list = re.split(r"\s+", esl_output_list[1])

        return float(lambda_list[2])