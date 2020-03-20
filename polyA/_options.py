import math

from argparse import ArgumentParser, Namespace
from typing import List, Optional, TextIO
from ._exceptions import ValidationException
from .constants import (
    DEFAULT_CHANGE_PROB,
    DEFAULT_CHUNK_SIZE,
    DEFAULT_SAME_PROB,
)
from .support_matrix import SupportMatrix, deserialize_support_matrix

# TODO: Can we use reflection to automate calling the _parse methods?
# TODO: Option for normal output target?
# TODO: Options for output formats?
# TODO: What if we want to handle different output matrices differently?
# TODO: Columns should be a list of columns to exclude for clarity
# TODO: Read options from a file by default or with a different script?
# TODO: Accept input from a pipe
# TODO: Add gap-ext and gap-init options and use defaults in constants


class Options:
    benchmark: bool
    change_prob: float
    same_prob: float
    chunk_size: int
    columns: Optional[List[int]]
    log_file: Optional[TextIO]
    logged_change_prob: float
    logged_same_prob: float
    support_matrix: SupportMatrix

    def __init__(self, args: Optional[List[str]] = None) -> None:
        if args is None:
            return

        parser = ArgumentParser(
            description="polyA adjudication tool", prog=__package__
        )
        parser.add_argument(
            "--benchmark",
            action="store_true",
            default=False,
            help="Collect and display benchmarking data",
        )
        parser.add_argument(
            "--change-prob",
            type=float,
            default=DEFAULT_CHANGE_PROB,
            help="Base probability of changing sequences",
        )
        parser.add_argument(
            "--chunk-size",
            type=int,
            default=DEFAULT_CHUNK_SIZE,
            help="Size of the window in base pairs analyzed together",
        )
        parser.add_argument(
            "--columns", type=str, help="Path list of column indices to run on",
        )
        parser.add_argument(
            "--log-file",
            type=str,
            help="File to log to or stdout, stderr, or none, default is stderr",
        )
        parser.add_argument(
            "--support", type=str, help="Path to a serialized support matrix",
        )

        namespace = parser.parse_args(args=args)
        self._parse_benchmark(namespace)
        self._parse_change_prob(namespace)
        self._parse_chunk_size(namespace)
        self._parse_columns(namespace)
        self._parse_log_file(namespace)
        self._parse_support_matrix(namespace)

    def _parse_benchmark(self, namespace: Namespace) -> None:
        self.benchmark = namespace.benchmark

    def _parse_change_prob(self, namespace: Namespace) -> None:
        self.change_prob = namespace.change_prob
        self.same_prob = 1.0 - self.change_prob

        if self.change_prob < 0 or self.change_prob > 1:
            raise ValidationException(
                f"--change-prob must be a probability in [0, 1], {self.change_prob} invalid"
            )

        self.logged_change_prob = math.log(self.change_prob)
        self.logged_same_prob = math.log(self.same_prob)

    def _parse_chunk_size(self, namespace: Namespace) -> None:
        self.chunk_size = namespace.chunk_size

    def _parse_columns(self, namespace: Namespace) -> None:
        if namespace.columns is None:
            self.columns = None
            return

        with open(namespace.columns) as columnsFile:
            columnsLines = columnsFile.readlines()
        self.columns = [int(c) for c in columnsLines]

    def _parse_log_file(self, namespace: Namespace) -> None:
        if namespace.log_file is None or namespace.log_file == "stderr":
            import sys

            self.log_file = sys.stderr
            return

        if namespace.log == "stdout":
            import sys

            self.log_file = sys.stdout
            return

        if namespace.log == "none":
            self.log_file = None
            return

        self.log_file = open(namespace.log, "a")

    def _parse_support_matrix(self, namespace: Namespace) -> None:
        if namespace.support is None:
            self.support_matrix = {}
            return

        with open(namespace.support) as supportFile:
            supportLines = supportFile.readlines()
        self.support_matrix = deserialize_support_matrix(supportLines)
