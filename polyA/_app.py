import logging
from sys import argv, stderr
from typing import List, IO

from ._options import Options
from ._runners import run_confidence, run_full
from .lambda_provider import EaselLambdaProvider
from .load_alignments import (
    load_alignments,
    load_alignment_tool,
    shard_overlapping_alignments,
)
from .output import Output
from .printers import Printer
from .prior_counts import read_prior_counts
from .substitution_matrix import load_substitution_matrices
from .ultra_provider import ApplicationUltraProvider, TandemRepeat


class AppError(RuntimeError):
    pass


def _configure_logging(opts: Options) -> None:
    log_file: IO[str]
    if opts.log_file_path != "":
        log_file = open(opts.log_file_path, "wt")
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

    # ----------------
    # File Conversions
    # ----------------

    if opts.cm_to_stockholm:
        from .converters.cm_to_stockholm import convert

        convert(opts.cm_to_stockholm, opts.stockholm_path, opts.matrix_path)

    if opts.rm_to_stockholm:
        from .converters.rm_to_stockholm import convert

        convert(opts.rm_to_stockholm, opts.stockholm_path, opts.matrix_path)

    if opts.cm_to_stockholm or opts.rm_to_stockholm:
        if not (opts.alignments_file_path and opts.sub_matrices_path):
            exit(0)

    # ----------------------------
    # Tandem Repeat initialization
    # ----------------------------

    # This will come back as an empty list if we don't need to worry about TRs.
    tandem_repeats = _configure_tandem_repeats(opts)

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

    with open(opts.alignments_file_path) as _align_tool_infile:
        alignment_tool: str = load_alignment_tool(_align_tool_infile)

    # Note: alignments from blast or HMMER are not set-up to use complexity
    # adjusted scoring
    if opts.complexity_adjustment and alignment_tool not in [
        "cross_match",
        "RepeatMasker",
    ]:
        raise AppError(
            f"cannot use complexity adjusted scoring with {alignment_tool}"
        )

    _lambda_provider = EaselLambdaProvider(opts.easel_path)
    with open(opts.sub_matrices_path) as _sub_matrices_file:
        sub_matrices = load_substitution_matrices(
            _sub_matrices_file, _lambda_provider, opts.complexity_adjustment
        )

    # -------------------------------------------------
    # Flags and parameters related to secondary outputs
    # -------------------------------------------------

    outputter = Output(opts.output_path, opts.output_to_file)

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

    # ---------------------
    # Set up general output
    # ---------------------

    results_file = outputter.get_results()
    printer = Printer(
        output_file=results_file,
        print_id=opts.ids,
        use_matrix_position=opts.matrix_position,
        use_sequence_position=opts.sequence_position,
    )
    printer.print_results_header()

    # ----------------------------------------------------------------
    # Loop through the alignment shards and process each independently
    # ----------------------------------------------------------------
    if opts.trans_penalty < 0:
        raise AppError(
            f"cannot use a transition penalty of {opts.trans_penalty}"
        )

    # FIXME: if shard gap is infinite, use all TRs (skip breaking them up)
    tr_start: int = 0
    _prev_start: int = -1
    _prev_stop: int = -1

    for index, chunk in enumerate(
            shard_overlapping_alignments(alignments, shard_gap=opts.shard_gap)
    ):
        chunk_start = chunk.start
        chunk_stop = chunk.stop
        # get TRs between chunk stop and start
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

        soda_viz_file = (
            outputter.get_soda(index) if opts.soda else None
        )
        printer.set_soda_files(soda_viz_file)

        (_last_start, _last_stop) = run_full(
            chunk.alignments,
            tandem_repeats_chunk,
            opts.chunk_size,
            opts.trans_penalty,
            sub_matrices,
            subfam_counts,
            chunk_start,
            chunk_stop,
            _prev_start,
            _prev_stop,
            printer,
        )
        _prev_start, _prev_stop = (
            _last_start,
            _last_stop,
        )

        if soda_viz_file is not None:
            soda_viz_file.close()

        # if soda_conf_file is not None:
        #     soda_conf_file.close()
