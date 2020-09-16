from logging import Formatter, INFO, Logger, StreamHandler, getLogger
from typing import List

from ._options import Options


def _configure_logging(options: Options) -> None:
    if options.log_file is not None:
        formatter = Formatter(fmt="{levelname} ({name}) - {message}", style="{")

        handler = StreamHandler(options.log_file)
        handler.setFormatter(formatter)

        logger = getLogger(__package__)
        logger.addHandler(handler)

        if options.benchmark:
            logger.setLevel(INFO)


def run(args: List[str]) -> None:
    options = Options(args)
    _configure_logging(options)
    print("foolish")

    # For now we just compute a probability matrix from a support
    # matrix because that's all we've ported to Python so far.

    # probMatrix, originMatrix = fill_prob_matrix(
    #     options.support_matrix,
    #     columns=options.columns,
    #     benchmark=options.benchmark,
    #     change_prob=options.change_prob,
    #     chunk_size=options.chunk_size,
    # )
    # serialize_prob_matrix(probMatrix, originMatrix)
