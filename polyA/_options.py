import math

from argparse import ArgumentParser, Namespace
from typing import List, Optional, TextIO
from ._exceptions import ValidationException
from .constants import (
    DEFAULT_CHANGE_PROB,
    DEFAULT_CHUNK_SIZE,
    DEFAULT_GAP_EXT,
    DEFAULT_GAP_INIT,
    DEFAULT_LAMBDA,
)
from .support_matrix import SupportMatrix, deserialize_support_matrix


# TODO: Can we use reflection to automate calling the _parse methods?
# TODO: Option for normal output target?
# TODO: Options for output formats?
# TODO: What if we want to handle different output matrices differently?
# TODO: Columns should be a list of columns to exclude for clarity
# TODO: Read options from a file by default or with a different script?
# TODO: Accept input from a pipe


class Options:
    """
    A typed container to hold program options and parameters.

    TODO: Decide how to document each option canonically

    >>> import sys
    >>> o = Options()
    >>> o.gap_ext == DEFAULT_GAP_EXT
    True
    >>> o.gap_init == DEFAULT_GAP_INIT
    True
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
            description="polyA adjudication tool", prog=__package__,
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
            "--gap-ext",
            type=int,
            help="TODO: Kaitlin",
            default=DEFAULT_GAP_EXT,
        )
        parser.add_argument(
            "--gap-init",
            type=int,
            help="TODO: Kaitlin",
            default=DEFAULT_GAP_INIT,
        )
        parser.add_argument(
            "--lambda",
            type=float,
            help="TODO: Kaitlin",
            default=DEFAULT_LAMBDA,
        )
        parser.add_argument(
            "--log-file",
            type=str,
            help="File to log to, or stdout, stderr, or none",
            default="stderr",
        )
        parser.add_argument(
            "--support", type=str, help="Path to a serialized support matrix",
        )

        namespace: Namespace
        if args is None:
            namespace = parser.parse_args(args=[])
        else:
            namespace = parser.parse_args(args=args)

        self._parse_all(namespace)

    def _parse_all(self, namespace: Namespace) -> None:
        self._parse_benchmark(namespace)
        self._parse_change_prob(namespace)
        self._parse_chunk_size(namespace)
        self._parse_columns(namespace)
        self._parse_gap_ext(namespace)
        self._parse_gap_init(namespace)
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
            columns_lines = columnsFile.readlines()
        self.columns = [int(c) for c in columns_lines]

    def _parse_gap_ext(self, namespace: Namespace) -> None:
        if namespace.gap_ext is None:
            self.gap_ext = DEFAULT_GAP_EXT
            return

        self.gap_ext = int(namespace.gap_ext)

    def _parse_gap_init(self, namespace: Namespace) -> None:
        if namespace.gap_init is None:
            self.gap_init = DEFAULT_GAP_INIT
            return

        self.gap_init = int(namespace.gap_init)

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

    def _parse_support_matrix(self, namespace: Namespace) -> None:
        if namespace.support is None:
            self.support_matrix = {}
            return

        with open(namespace.support) as supportFile:
            support_lines = supportFile.readlines()
        self.support_matrix = deserialize_support_matrix(support_lines)
