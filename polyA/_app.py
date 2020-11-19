import re
from getopt import getopt
from sys import argv
from typing import Dict, List

from polyA._runners import run_confidence, run_full
from polyA import *

helpMessage: str = f"""
usage: {argv[0]} alignFile subMatrixFile\n
ARGUMENTS
    --gap-init [-25]
    --gap-ext [-5]
    --lambda [will calculate from substitution matrix if not included]
    --chunk-size (must be odd) [31]
    --esl-path <path to Easel>
    --confidence - output confidence for a single annoation without running whole algorithm
    --prior-counts <path to file>
    --ultra-path <path to ultra>
    --seq-file <path to file>
    --ultra-output <path to file>
    --viz <output path> - prints output format for SODA visualization
    --heatmap <output path> - prints probability file for input into heatmap

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
            "viz=",
            "heatmap=",
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
    lambdaa = float(opts["--lambda"]) if "--lambda" in opts else 0.0
    esl_path = str(opts["--esl-path"]) if "--esl-path" in opts else ""
    chunk_size = (
        int(opts["--chunk-size"])
        if "--chunk-size" in opts
        else DEFAULT_CHUNK_SIZE
    )

    infile_prior_counts = (
        open(opts["--prior-counts"], "r") if "--prior-counts" in opts else None
    )

    outfile_viz_path = str(opts["--viz"]) if "--viz" in opts else ""
    outfile_viz = open(outfile_viz_path, "w") if outfile_viz_path else None
    outfile_conf = (
        open(outfile_viz_path + ".json", "w") if outfile_viz_path else None
    )

    outfile_heatmap_path = str(opts["--heatmap"]) if "--heatmap" in opts else ""
    outfile_heatmap = (
        open(outfile_heatmap_path, "w") if outfile_heatmap_path else None
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

    # input is alignment file of hits region and substitution matrix
    infile: str = args[0]
    infile_matrix: str = args[1]

    # Other open was moved down to where we load the alignments file
    with open(infile_matrix) as _infile_matrix:
        in_matrix: List[str] = _infile_matrix.readlines()

    # if lambda isn't included at command line, run esl_scorematrix to calculate it from scorematrix
    if not lambdaa:
        provider = EaselLambdaProvider(esl_path, infile_matrix)
        lambdaa = provider()

    # reads in the score matrix from file and stores in dict that maps 'char1char2' to the score from the
    # input substitution matrix - ex: 'AA' = 8
    sub_matrix: Dict[str, int] = {}

    # add all ambiguity codes just incase matrix doesnt have them
    nucleotide_codes = "AGCTYRWSKMDVHBXN."
    for code in nucleotide_codes:
        for code2 in nucleotide_codes:
            sub_matrix[code + code2] = 0

    line = in_matrix[0]
    line = re.sub(r"^\s+", "", line)
    line = re.sub(r"\s+$", "", line)
    chars = re.split(r"\s+", line)

    count: int = 0
    for line in in_matrix[1:]:
        line = re.sub(r"^\s+", "", line)
        line = re.sub(r"\s+$", "", line)
        sub_scores = re.split(r"\s+", line)
        for i in range(len(sub_scores)):
            sub_matrix[chars[count] + chars[i]] = int(sub_scores[i])
        count += 1
    sub_matrix[".."] = 0

    subfam_counts = (
        read_prior_counts(infile_prior_counts, using_tr)
        if infile_prior_counts is not None
        else {}
    )

    with open(infile) as _infile:
        alignments = list(load_alignments(_infile))

    if confidence_flag:
        run_confidence(alignments, lambdaa=lambdaa)
        exit()

    run_full(
        alignments,
        tandem_repeats,
        chunk_size,
        gap_ext,
        gap_init,
        lambdaa,
        outfile_conf,
        outfile_viz,
        outfile_heatmap,
        print_matrix_pos,
        print_seq_pos,
        sub_matrix,
        subfam_counts,
    )

    for fp in [
        infile_prior_counts,
        outfile_conf,
        outfile_viz,
        outfile_heatmap,
    ]:
        if fp is not None:
            fp.close()
