from argparse import ArgumentParser, Namespace
from typing import List, Optional, TextIO
from .support_matrix import SupportMatrix, deserialize_support_matrix

# TODO: Can we use reflection to automate calling the _parse methods?
# TODO: Option for normal output target?
# TODO: Options for output formats?
# TODO: What if we want to handle different output matrices differently?
# TODO: Columns should be a list of columns to exclude for clarity
# TODO: Read options from a file by default or with a different script?
# TODO: Accept input from a pipe
class Options:
    benchmark: bool
    columns: Optional[List[int]]
    log_target: Optional[TextIO]
    support_matrix: SupportMatrix

    def __init__(self, args: Optional[List[str]] = None) -> None:
        if args is None:
            return

        parser = ArgumentParser(description="polyA adjudication tool")
        parser.add_argument(
            "--benchmark",
            action="store_true",
            default=False,
            help="Collect and display benchmarking data",
        )
        parser.add_argument(
            "--columns", type=str, help="Path list of column indices to run on"
        )
        parser.add_argument(
            "--log",
            type=str,
            help="File to log to or stdout, stderr, or none, default is stderr",
        )
        parser.add_argument(
            "--support", type=str, help="Path to a serialized support matrix"
        )

        namespace = parser.parse_args(args=args)
        self._parse_benchmark(namespace)
        self._parse_columns(namespace)
        self._parse_log(namespace)
        self._parse_support_matrix(namespace)

    def _parse_benchmark(self, namespace: Namespace) -> None:
        self.benchmark = namespace.benchmark

    def _parse_columns(self, namespace: Namespace) -> None:
        if namespace.columns is None:
            self.columns = None
            return

        with open(namespace.columns) as columnsFile:
            columnsLines = columnsFile.readlines()
        self.columns = [int(c) for c in columnsLines]

    def _parse_log(self, namespace: Namespace) -> None:
        if namespace.log is None or namespace.log == "stderr":
            import sys

            self.log_target = sys.stderr
            return

        if namespace.log == "stdout":
            import sys

            self.log_target = sys.stdout
            return

        if namespace.log == "none":
            self.log_target = None
            return

        self.log_target = open(namespace.log, "a")

    def _parse_support_matrix(self, namespace: Namespace) -> None:
        if namespace.support is None:
            self.support_matrix = {}
            return

        with open(namespace.support) as supportFile:
            supportLines = supportFile.readlines()
        self.support_matrix = deserialize_support_matrix(supportLines)
