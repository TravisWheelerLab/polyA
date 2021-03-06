import logging
from sys import argv, stderr, stdout
from typing import List

from ._options import Options
from ._runners import run_confidence, run_full
from .lambda_provider import EaselLambdaProvider
from .load_alignments import (
    load_alignments,
    shard_overlapping_alignments,
)
from .output import Output
from .prior_counts import read_prior_counts
from .substitution_matrix import load_substitution_matrices
from .ultra_provider import ApplicationUltraProvider, TandemRepeat


class AppError(RuntimeError):
    pass


def _configure_logging(opts: Options) -> None:
    if opts.log_file_path != "":
        log_file = open(opts.log_file_path, "w")
    else:
        log_file = stderr
    handler = logging.StreamHandler(log_file)

    formatter = logging.Formatter(
        fmt="{levelname} ({name}) - {message}", style="{"
    )
    handler.setFormatter(formatter)

    logger = logging.getLogger(__package__)
    logger.addHandler(handler)

    if opts.log_level == "debug":
        logger.setLevel(logging.DEBUG)
    elif opts.log_level == "verbose":
        logger.setLevel(logging.INFO)
    elif opts.log_level == "normal":
        logger.setLevel(logging.WARNING)
    elif opts.log_level == "quiet":
        logger.setLevel(logging.CRITICAL)


def _configure_tandem_repeats(opts: Options) -> List[TandemRepeat]:
    if opts.ultra_data_path and opts.sequence_file_path:
        # This is ambiguous so we bail out
        raise AppError("cannot specify both ultra data and sequence files")

    tandem_repeats = []

    if opts.ultra_data_path or opts.sequence_file_path:
        # Run ULTRA or load TRs from an output file
        provider = ApplicationUltraProvider(
            opts.sequence_file_path,
            opts.ultra_data_path,
            opts.ultra_path,
        )
        ultra_output = provider()
        tandem_repeats = ultra_output.tandem_repeats

    return tandem_repeats


def run():
    opts = Options(argv[1:])
    _configure_logging(opts)

    # ----------------------------
    # Tandem Repeat initialization
    # ----------------------------

    tandem_repeats = _configure_tandem_repeats(
        opts
    )  # either has stuff or doesn't

    # -----------------
    # Sub-family counts
    # -----------------

    infile_prior_counts = (
        open(opts.prior_counts_path, "r") if opts.prior_counts_path else None
    )

    subfam_counts = (
        read_prior_counts(infile_prior_counts, bool(tandem_repeats))
        if infile_prior_counts is not None
        else {}
    )

    # ----------------------------
    # Load the substitution matrix
    # ----------------------------

    _lambda_provider = EaselLambdaProvider(opts.easel_path)
    with open(opts.sub_matrices_path) as _sub_matrices_file:
        sub_matrices = load_substitution_matrices(
            _sub_matrices_file, _lambda_provider
        )

    # -------------------------------------------------
    # Flags and parameters related to secondary outputs
    # -------------------------------------------------

    outputter = Output(opts.output_path)

    # -----------------------------
    # Load alignments to operate on
    # -----------------------------

    with open(opts.alignments_file_path) as _infile:
        alignments = list(load_alignments(_infile))

    # --------------------------
    # Run confidence calculation
    # --------------------------

    lambda_values = [sub_matrices[a.sub_matrix_name].lamb for a in alignments]

    if opts.confidence:
        run_confidence(
            alignments,
            lambs=lambda_values,
        )
        exit()

    # ----------------------------------------------------------------
    # Loop through the alignment shards and process each independently
    # ----------------------------------------------------------------

    # FIXME: if shard gap is infinite, use all TRs (skip breaking them up)
    stdout.write("start\tstop\tID\tname\n")
    stdout.write("----------------------------------------\n")
    tr_start: int = 0
    tr_end: int = 0
    _prev_start: int = -1
    _prev_stop: int = -1
    for index, chunk in enumerate(
        shard_overlapping_alignments(alignments, shard_gap=opts.shard_gap)
    ):
        chunk_start = chunk.start
        chunk_stop = chunk.stop
        # get TRs between chunk stop and start
        tandem_repeats_chunk: List[TandemRepeat] = []
        tr_end = tr_start
        while tr_end < len(tandem_repeats):
            tr = tandem_repeats[tr_end]
            if tr.start <= chunk_stop:
                tr_end += 1
            else:
                break
        tandem_repeats_chunk = tandem_repeats[tr_start:tr_end]
        tr_start = tr_end
        if tr_start > 0 and tandem_repeats[tr_start - 1].stop > chunk_stop:
            tr_start -= 1

        soda_viz_file, soda_conf_file = (
            outputter.get_soda(index) if opts.soda else (None, None)
        )

        (_last_start, _last_stop) = run_full(
            chunk.alignments,
            tandem_repeats_chunk,
            opts.chunk_size,
            soda_viz_file,
            soda_conf_file,
            opts.matrix_position,
            opts.sequence_position,
            sub_matrices,
            subfam_counts,
            chunk_start,
            chunk_stop,
            _prev_start,
            _prev_stop,
        )
        _prev_start, _prev_stop = (
            _last_start,
            _last_stop,
        )

        if soda_viz_file is not None:
            soda_viz_file.close()
        if soda_conf_file is not None:
            soda_conf_file.close()
