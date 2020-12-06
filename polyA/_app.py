from getopt import getopt
import logging
from sys import argv

from polyA._runners import run_confidence, run_full
from polyA import *

logging.root.addHandler(logging.StreamHandler())
logging.root.setLevel(logging.DEBUG)

helpMessage: str = f"""
usage: {argv[0]} alignFile subMatrixFile\n
ARGUMENTS
    --gap-init [-25]
    --gap-ext [-5]
    --chunk-size (must be odd) [31]
    --esl-path <path to Easel>
    --confidence - output confidence for a single annotation without running whole algorithm

    --output-path - path to a directory where output files will be written
    --heatmap - output a heatmap file
    --soda - output files to produce a SODA visualization

    --prior-counts <path to file>
    --ultra-path <path to ultra>
    --seq-file <path to file>
    --ultra-output <path to file>
    --shard-gap <int> - allowed gap between sequences in the same shard

OPTIONS
    --help - display help message
    --matrix-pos - prints output in terms of matrix position
    --seq-pos - prints output in terms of target sequence position
"""


def run():
    raw_opts, args = getopt(
        argv[1:],
        "",
        [
            "gap-init=",
            "gap-ext=",
            "skipScore=",
            "lambda=",
            "esl-path=",
            "confidence",
            "chunk-size=",
            "prior-counts=",
            "output-path=",
            "soda",
            "heatmap",
            "ultra-path=",
            "seq-file=",
            "ultra-output=",
            "help",
            "matrix-pos",
            "seq-pos",
        ],
    )
    opts = dict(raw_opts)

    help_flag = "--help" in opts
    if help_flag:
        print(helpMessage)
        exit(0)

    gap_init = (
        int(opts["--gap-init"]) if "--gap-init" in opts else DEFAULT_GAP_INIT
    )
    gap_ext = int(opts["--gap-ext"]) if "--gap-ext" in opts else DEFAULT_GAP_EXT
    esl_path = str(opts["--esl-path"]) if "--esl-path" in opts else ""
    chunk_size = (
        int(opts["--chunk-size"])
        if "--chunk-size" in opts
        else DEFAULT_CHUNK_SIZE
    )

    infile_prior_counts = (
        open(opts["--prior-counts"], "r") if "--prior-counts" in opts else None
    )

    # Run ULTRA or load TRs from an output file
    tandem_repeats = []
    # TODO: Assume ULTRA is in the PATH if not given
    # Right now we assume that if the user provided either an ultra
    # path OR a path to an ultra output file then they want to run
    # the TR code. In the future we will want a better way to indicate
    # that since we will assume ultra is in the PATH.
    ultra_path = str(opts["--ultra-path"]) if "--ultra-path" in opts else ""
    ultra_output_path = (
        str(opts["--ultra-output"]) if "--ultra-output" in opts else ""
    )
    using_tr = bool(ultra_output_path) or bool(ultra_path)
    if using_tr:
        seq_file_path = str(opts["--seq-file"]) if "--seq-file" in opts else ""
        provider = ApplicationUltraProvider(
            seq_file_path, ultra_output_path, ultra_path
        )
        ultra_output = provider()
        tandem_repeats = ultra_output.tandem_repeats

    print_matrix_pos = "--matrix-pos" in opts
    print_seq_pos = "--seq-pos" in opts
    confidence_flag = "--confidence" in opts

    shard_gap = (
        int(opts["--shard-gap"]) if "--shard-gap" in opts else DEFAULT_SHARD_GAP
    )

    # input is alignment file of hits region and substitution matrix
    infile: str = args[0]

    subfam_counts = (
        read_prior_counts(infile_prior_counts, using_tr)
        if infile_prior_counts is not None
        else {}
    )

    # ----------------------------
    # Load the substitution matrix
    # ----------------------------

    lambda_provider = EaselLambdaProvider(esl_path)

    sub_matrix_path: str = args[1]
    with open(sub_matrix_path) as _sub_matrix_file:
        sub_matrices = load_substitution_matrices(
            _sub_matrix_file, lambda_provider
        )

    # -------------------------------------------------
    # Flags and parameters related to secondary outputs
    # -------------------------------------------------

    soda_flag = "--soda" in opts
    heatmap_flag = "--heatmap" in opts

    output_path = (
        opts["--output-path"] if "--output-path" in opts else "polya-output"
    )
    outputter = Output(output_path)

    # -----------------------------
    # Load alignments to operate on
    # -----------------------------

    with open(infile) as _infile:
        alignments = list(load_alignments(_infile))

    lambda_values = [sub_matrices[a.sub_matrix_name].lamb for a in alignments]

    # --------------------------
    # Run confidence calculation
    # --------------------------

    if confidence_flag:
        run_confidence(
            alignments,
            lambs=lambda_values,
        )
        exit()

    # ----------------------------------------------------------------
    # Loop through the alignment shards and process each independently
    # ----------------------------------------------------------------

    for index, chunk in enumerate(
        shard_overlapping_alignments(alignments, shard_gap=shard_gap)
    ):
        soda_viz_file, soda_conf_file = (
            outputter.get_soda(index) if soda_flag else (None, None)
        )
        heatmap_file = outputter.get_heatmap(index) if heatmap_flag else None

        run_full(
            chunk,
            tandem_repeats,
            chunk_size,
            gap_ext,
            gap_init,
            lambda_values,
            soda_viz_file,
            soda_conf_file,
            heatmap_file,
            print_matrix_pos,
            print_seq_pos,
            sub_matrices,
            subfam_counts,
        )

        if soda_viz_file is not None:
            soda_viz_file.close()
        if soda_conf_file is not None:
            soda_conf_file.close()
        if heatmap_file is not None:
            heatmap_file.close()
