import math
from argparse import ArgumentParser, Namespace
from typing import List, Optional, TextIO

from _exceptions import ValidationException
from alignment import Alignment
from constants import (
    DEFAULT_CHANGE_PROB,
    DEFAULT_CHUNK_SIZE,
    DEFAULT_GAP_EXTEND,
    DEFAULT_GAP_START,
    DEFAULT_LAMBDA,
)
from edges import edges
from load_alignments import load_alignments
from pad_sequences import pad_sequences
from polyA import SubstitutionMatrix
from support_matrix import SupportMatrix, deserialize_support_matrix


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
    >>> o.gap_ext == DEFAULT_GAP_EXTEND
    True
    >>> o.gap_init == DEFAULT_GAP_START
    True
    >>> o.log_file.name == sys.stderr.name
    True
    >>> o = Options(['--log-file', 'foo.txt'])
    >>> o.log_file.name
    'foo.txt'
    """

    alignments: List[Alignment]
    benchmark: bool
    change_prob: float
    chunk_size: int
    columns: Optional[List[int]]
    edge_start: int
    edge_stop: int
    gap_ext: int
    gap_init: int
    lambda: float
    log_file: Optional[TextIO]
    logged_change_prob: float
    logged_skip_change_prob: float
    logged_skip_same_prob: float
    same_prob: float
    # TODO: Write a _parse method that loads this
    substitution_matrix: SubstitutionMatrix
    support_matrix: SupportMatrix

    def __init__(self, args: Optional[List[str]] = None) -> None:
        parser = ArgumentParser(
            description="polyA adjudication tool", prog=__package__,
        )

        parser.add_argument(
            "alignments-path", help="Path to the alignments file to process",
        )
        parser.add_argument(
            "substitution-path", help="Path to the substitution matrix",
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
            default=DEFAULT_GAP_EXTEND,
        )
        parser.add_argument(
            "--gap-init",
            type=int,
            help="TODO: Kaitlin",
            default=DEFAULT_GAP_START,
        )
        parser.add_argument(
            "--lambda-value", type=float, help="TODO: Kaitlin", default=None,
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
        self._parse_alignments_path(namespace)
        self._parse_benchmark(namespace)
        self._parse_change_prob(namespace)
        self._parse_chunk_size(namespace)
        self._parse_columns(namespace)
        self._parse_gap_ext(namespace)
        self._parse_gap_init(namespace)
        self._parse_log_file(namespace)
        self._parse_support_matrix(namespace)

    def _parse_alignments_path(self, namespace: Namespace) -> None:
        path: str = namespace.alignments_path
        # TODO: Check for stdin / pipe
        with open(path, "r") as file:
            alignments = load_alignments(file)

        # TODO: Should this go here or in load_alignments?
        pad_sequences(alignments)

        edge_start, edge_stop = edges(alignments)
        self.edge_start = edge_start
        self.edge_stop = edge_stop

        self.alignments = alignments

    def _parse_benchmark(self, namespace: Namespace) -> None:
        self.benchmark = namespace.benchmark

    def _parse_change_prob(self, namespace: Namespace) -> None:
        # Note that _parse_alignments_path must have already run
        # at this point because we use the alignments

        self.change_prob = namespace.change_prob
        self.same_prob = 1.0 - self.change_prob

        if self.change_prob < 0 or self.change_prob > 1:
            raise ValidationException(
                f"--change-prob must be a probability in [0, 1], {self.change_prob} invalid"
            )

        num_alignments = len(self.alignments)
        self.logged_change_prob = math.log(
            self.change_prob / (num_alignments - 1)
        )
        self.logged_same_prob = math.log(self.same_prob)

        self.logged_skip_same_prob = self.logged_change_prob / 10.0
        self.logged_skip_change_prob = self.logged_change_prob / 2.0

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
            self.gap_ext = DEFAULT_GAP_EXTEND
            return

        self.gap_ext = int(namespace.gap_ext)

    def _parse_gap_init(self, namespace: Namespace) -> None:
        if namespace.gap_init is None:
            self.gap_init = DEFAULT_GAP_START
            return

        self.gap_init = int(namespace.gap_init)

    def _parse_lambda(self, namespace: Namespace) -> None:
        if namespace.lambda_value is None:
            self.lambda_value = namespace.lambda_value
            # TODO: Add option to run esl_scorematrix
            return

        self.lambda_value = float(namespace.lambda_value)

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
