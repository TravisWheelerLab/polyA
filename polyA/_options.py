from argparse import ArgumentParser, Namespace
from typing import List, Literal, Optional

from . import __version__
from .constants import (
    DEFAULT_CHUNK_SIZE,
    DEFAULT_SHARD_GAP,
)


class Options:
    """
    A typed container to hold program options and parameters.

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

    log_file_path: str
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

        parser.add_argument(
            "alignments_file_path",
            metavar="ALIGNMENTS",
            help="alignments file in Stockholm format",
        )
        parser.add_argument(
            "sub_matrices_path",
            metavar="MATRICES",
            help="substitution matrices file in PolyA matrix format",
        )

        parser.add_argument(
            "-v",
            "--version",
            action="version",
            version=__version__,
            help="show version and exit",
        )

        parser.add_argument(
            "--chunk-size",
            type=int,
            default=DEFAULT_CHUNK_SIZE,
            help="size of the window in base pairs analyzed together",
        )
        parser.add_argument(
            "--confidence",
            action="store_true",
            default=False,
            help="run the confidence calculation and then exit",
        )
        parser.add_argument(
            "--prior-counts",
            metavar="FILE",
            default="",
            help="file containing query genomic counts",
        )
        parser.add_argument(
            "--shard-gap",
            type=int,
            default=DEFAULT_SHARD_GAP,
            help="maximum alignment gap before sharding occurs",
        )
        parser.add_argument(
            "--sequences",
            metavar="SEQS",
            default="",
            help="fasta file for running ULTRA",
        )
        parser.add_argument(
            "--ultra-data",
            metavar="FILE",
            default="",
            help="file of the output from ULTRA",
        )

        parser.add_argument(
            "--easel-path",
            metavar="BIN",
            default="",
            help="path to the esl_scorematrix program, if necessary (assumed to be in PATH)",
        )
        parser.add_argument(
            "--ultra-path",
            metavar="BIN",
            default="ultra",
            help="path to the ULTRA binary to use, if necessary (assumed to be in PATH)",
        )

        parser.add_argument(
            "--log-file",
            metavar="LOG",
            default="",
            help="file to store log output in, defaults to stderr",
        )
        parser.add_argument(
            "--log-level",
            metavar="LEVEL",
            choices=["debug", "verbose", "normal", "quiet"],
            help="logging level to use, 'debug' is the most noisy",
        )
        parser.add_argument(
            "--matrix-position",
            action="store_true",
            default=False,
            help="produce output in terms of the matrix position",
        )
        parser.add_argument(
            "--output-path",
            metavar="PATH",
            default=".",
            help="directory to write output files to, defaults to working directory",
        )
        parser.add_argument(
            "--sequence-position",
            action="store_true",
            default=False,
            help="produce output in terms of the target sequence position",
        )
        parser.add_argument(
            "--soda",
            action="store_true",
            default=False,
            help="write a SODA visualization file to the output directory",
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

        self.log_file_path = namespace.log_file
        self.log_level = namespace.log_level
        self.matrix_position = namespace.matrix_position
        self.output_path = namespace.output_path
        self.sequence_position = namespace.sequence_position
        self.soda = namespace.soda
