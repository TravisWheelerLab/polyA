from argparse import ArgumentParser, Namespace
from typing import List, Literal, Optional

from .constants import (
    DEFAULT_CHUNK_SIZE,
    DEFAULT_SHARD_GAP,
)


class Options:
    """
    A typed container to hold program options and parameters.

    TODO(George): Option for normal output target?
    TODO(George): Options for output formats?
    TODO(George): Read options from a file by default or with a different script?
    TODO(George): Accept input from a pipe?

    >>> o = Options()
    >>> o.log_file_path
    ''
    >>> o = Options(["", "", "--log-file", "foo.txt"])
    >>> o.log_file_path
    'foo.txt'
    """

    alignments_file_path: str
    sub_matrices_path: str

    # -----------------------
    # Algorithm Configuration
    # -----------------------

    chunk_size: int
    confidence: bool
    prior_counts_path: str
    shard_gap: int
    sequence_file_path: str
    ultra_data_path: str

    # -------------------
    # Helper applications
    # -------------------

    easel_path: str
    ultra_path: str

    # --------------------
    # Output configuration
    # --------------------

    heatmap: bool
    log_file_path: str
    # TODO(George): Might as well make this the actual log level
    log_level: Literal["debug", "verbose", "normal", "quiet"]
    matrix_position: bool
    output_path: str
    sequence_position: bool
    soda: bool

    def __init__(self, args: Optional[List[str]] = None) -> None:
        parser = ArgumentParser(
            description="PolyA sequence adjudication tool",
            prog=__package__,
        )

        # FIXME(George): Make sure these each only consume one argument
        parser.add_argument(
            "alignments_file_path",
            metavar="FILE",
            help="Alignments file in Stockholm format",
        )
        parser.add_argument(
            "sub_matrices_path",
            metavar="FILE",
            help="Substitution matrices file in PolyA matrix format",
        )

        parser.add_argument(
            "--chunk-size",
            type=int,
            default=DEFAULT_CHUNK_SIZE,
            help="Size of the window in base pairs analyzed together",
        )
        parser.add_argument(
            "--confidence",
            action="store_true",
            default=False,
            help="Run the confidence calculation and then exit",
        )
        parser.add_argument(
            "--prior-counts",
            metavar="FILE",
            default="",
            help="TODO(Kaitlin)",
        )
        parser.add_argument(
            "--shard-gap",
            type=int,
            default=DEFAULT_SHARD_GAP,
            help="Maximum alignment gap before sharding occurs",
        )
        # TODO(Aubrey): Decide on a better name for this one, if possible
        parser.add_argument(
            "--sequences",
            metavar="FILE",
            default="",
            help="TODO(Aubrey)",
        )
        parser.add_argument(
            "--ultra-data",
            metavar="FILE",
            default="",
            help="TODO(Audrey)",
        )

        parser.add_argument(
            "--easel-path",
            metavar="BIN",
            default="esl_scorematrix",
            help="Path to the esl_scorematrix program, if necessary (assumed to be in PATH)",
        )
        parser.add_argument(
            "--ultra-path",
            metavar="BIN",
            default="ultra",
            help="Path to the ULTRA binary to use, if necessary (assumed to be in PATH)",
        )

        parser.add_argument(
            "--heatmap",
            action="store_true",
            default=False,
            help="Write a heatmap file to the output directory",
        )
        parser.add_argument(
            "--log-file",
            metavar="FILE",
            default="",
            help="File to store log output in, defaults to stderr",
        )
        parser.add_argument(
            "--log-level",
            metavar="LEVEL",
            choices=["debug", "verbose", "normal", "quiet"],
            help="Logging level to use, 'debug' is the most noisy",
        )
        parser.add_argument(
            "--matrix-position",
            action="store_true",
            default=False,
            help="Produce output in terms of the matrix position",
        )
        parser.add_argument(
            "--output-path",
            metavar="PATH",
            default=".",
            help="Directory to write output files to, defaults to working directory",
        )
        parser.add_argument(
            "--sequence-position",
            action="store_true",
            default=False,
            help="Produce output in terms of the target sequence position",
        )
        parser.add_argument(
            "--soda",
            action="store_true",
            default=False,
            help="Write a SODA visualization file to the output directory",
        )

        namespace: Namespace
        if args is None:
            namespace = parser.parse_args(args=["", ""])
        else:
            namespace = parser.parse_args(args=args)

        self.alignments_file_path = namespace.alignments_file_path
        self.sub_matrices_path = namespace.sub_matrices_path

        self.chunk_size = namespace.chunk_size
        self.confidence = namespace.confidence
        self.prior_counts_path = namespace.prior_counts
        self.shard_gap = namespace.shard_gap
        self.sequence_file_path = namespace.sequences
        self.ultra_data_path = namespace.ultra_data

        self.easel_path = namespace.easel_path
        self.ultra_path = namespace.ultra_path

        self.heatmap = namespace.heatmap
        self.log_file_path = namespace.log_file
        self.log_level = namespace.log_level
        self.matrix_position = namespace.matrix_position
        self.output_path = namespace.output_path
        self.sequence_position = namespace.sequence_position
        self.soda = namespace.soda
