import math

from argparse import ArgumentParser, Namespace
from typing import List, Optional, TextIO

from polyA._exceptions import ValidationException
from polyA.constants import (
    DEFAULT_CHUNK_SIZE,
)
from polyA.matrices import SupportMatrix


# TODO: George - Can we use reflection to automate calling the _parse methods?
# TODO: George - Option for normal output target?
# TODO: George - Options for output formats?
# TODO: George - What if we want to handle different output matrices differently?
# TODO: George - Columns should be a list of columns to exclude for clarity
# TODO: George - Read options from a file by default or with a different script?
# TODO: George - Accept input from a pipe


class Options:
    """
    A typed container to hold program options and parameters.

    TODO: George - Decide how to document each option canonically

    >>> import sys
    >>> o = Options()
    >>> o.log_file.name == sys.stderr.name
    True
    >>> o = Options(['--log-file', 'foo.txt'])
    >>> o.log_file.name
    'foo.txt'
    """

    benchmark: bool
    change_prob: float
    chunk_size: int
    columns: Optional[List[int]]
    gap_ext: int
    gap_init: int
    log_file: Optional[TextIO]
    logged_change_prob: float
    logged_same_prob: float
    same_prob: float
    support_matrix: SupportMatrix

    def __init__(self, args: Optional[List[str]] = None) -> None:
        parser = ArgumentParser(
            description="polyA adjudication tool",
            prog=__package__,
        )

        parser.add_argument(
            "--benchmark",
            action="store_true",
            default=False,
            help="Collect and display benchmarking data",
        )
        parser.add_argument(
            "--chunk-size",
            type=int,
            default=DEFAULT_CHUNK_SIZE,
            help="Size of the window in base pairs analyzed together",
        )
        parser.add_argument(
            "--columns",
            type=str,
            help="Path list of column indices to run on",
        )
        parser.add_argument(
            "--log-file",
            type=str,
            help="File to log to, or stdout, stderr, or none",
            default="stderr",
        )
        parser.add_argument(
            "--support",
            type=str,
            help="Path to a serialized support matrix",
        )

        namespace: Namespace
        if args is None:
            namespace = parser.parse_args(args=[])
        else:
            namespace = parser.parse_args(args=args)

        self._parse_all(namespace)

    def _parse_all(self, namespace: Namespace) -> None:
        self._parse_benchmark(namespace)
        self._parse_chunk_size(namespace)
        self._parse_columns(namespace)
        self._parse_log_file(namespace)

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
            columns_lines = columnsFile.readlines()
        self.columns = [int(c) for c in columns_lines]

    def _parse_log_file(self, namespace: Namespace) -> None:
        if namespace.log_file is None or namespace.log_file == "stderr":
            import sys

            self.log_file = sys.stderr
            return

        if namespace.log_file == "stdout":
            import sys

            self.log_file = sys.stdout
            return

        if namespace.log_file == "none":
            self.log_file = None
            return

        self.log_file = open(namespace.log_file, "a")
